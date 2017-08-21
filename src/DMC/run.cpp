#include "DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::run () {
  // initiate/reset FeynmanDiagram
  this->FD = FeynmanDiagram{this->initialExternalMomentum, this->dt, this->alpha, this->mu};

  // initiate/reset histogram
  unsigned long long int itrNum = 0;
  this->N0 = 0;
  this->hist = Array<unsigned long long int, Dynamic, Dynamic>::Zero(this->Np, this->Nt);

  // save an empty file
  if (this->worldRank == 0) {
    this->write2file(this->hist, this->N0, itrNum, 0);
  }

  // used for temporary saving result
  time_t startTime = time(NULL);

  // MC instead of DMC for S1 at small times
  // this value is overwritten
  this->MCvsDMCboundary = 0;
  if (this->minDiagramOrder <= 1 && ! this->externalLegs) {
    this->firstOrderSelfEnergyMC(this->numSecsDoingMCperCorePerBoldItr);
  }

  // specify the relative probability of choosing a specific update function
  multimap<unsigned int, void (DiagrammaticMonteCarlo::*)(double)> updateMethods;
  if (this->bold && this->boldIteration > 0) {
    updateMethods = {
      {1, &DiagrammaticMonteCarlo::BOLDshiftVertexPosition},
      {1, &DiagrammaticMonteCarlo::swapPhononConnections},
      {1, &DiagrammaticMonteCarlo::changeInternalPhononMomentumDirection},
      {1, &DiagrammaticMonteCarlo::changeInternalPhononMomentumMagnitude},
      {1, &DiagrammaticMonteCarlo::BOLDraiseOrder}, // <- These two must have the same probability
      {1, &DiagrammaticMonteCarlo::BOLDlowerOrder}, // <-
      {1, &DiagrammaticMonteCarlo::BOLDchangeDiagramLength},
      {1, &DiagrammaticMonteCarlo::BOLDchangeDiagramLengthComplex},
      {(this->fixedExternalMomentum ? 0 : 1), &DiagrammaticMonteCarlo::changeExternalMomentumMagnitude}
    };
  } else {
    updateMethods = {
      {1, &DiagrammaticMonteCarlo::shiftVertexPosition},
      {1, &DiagrammaticMonteCarlo::swapPhononConnections},
      {1, &DiagrammaticMonteCarlo::changeInternalPhononMomentumDirection},
      {1, &DiagrammaticMonteCarlo::changeInternalPhononMomentumMagnitude},
      {2, &DiagrammaticMonteCarlo::changeDiagramOrder},
      {1, &DiagrammaticMonteCarlo::changeDiagramLength},
      {1, &DiagrammaticMonteCarlo::changeDiagramLengthComplex},
      {(this->fixedExternalMomentum ? 0 : 1), &DiagrammaticMonteCarlo::changeExternalMomentumMagnitude}
    };
  }

  // vector which is going to contain the specified quantity of update functions
  vector<void (DiagrammaticMonteCarlo::*)(double)> chooseUpdateMethod;

  // populate vector
  for (auto updateMethod = updateMethods.begin(); updateMethod != updateMethods.end(); updateMethod++) {
    for (unsigned int i = 0; i != updateMethod->first; i++) {
      chooseUpdateMethod.push_back(updateMethod->second);
    }
  }

  // take into account time spent reaching the start configuration
  const clock_t beginTime = clock();
  
  // to reach a random start connfiguration
  for (unsigned long long int i = 0; i < this->untilStart; ++i) {
    // choose update operation on random
    auto updateMethod = chooseUpdateMethod[this->Uint(0, chooseUpdateMethod.size() - 1)];
    // update diagram
    (this->*updateMethod)(this->param);
  }

  // // number of iterations for each process
  // unsigned long long int localNumIterations = this->numIterations/this->worldSize;
  // long long int rest = this->numIterations - localNumIterations*this->worldSize;
  // if (this->worldRank < rest) localNumIterations++;

  // // when to do a temporary save
  // const unsigned long long int numItrBetweenTempSaves = localNumIterations/(this->numTempDMCsaves + 1);
  // unsigned long long int numItrAfterPrevTempSave = 0;

  double secsPerTempSave = numSecsDoingDMCperCorePerBoldItr/(this->numTempDMCsaves + 1);

  // a fixed number of iterations we do before looking at the time difference
  // this since looking at the time difference is a rather complex computation
  unsigned int numItrsPerTimesCheck = 1000000;

  // first save has already occured, thus initialize to 0
  unsigned int tempSaveNum = 0;

  // main loop
  // for (itrNum = 0; itrNum < localNumIterations; itrNum++) {


  unsigned long long int totMeanTimes = 0;
  double meanTime = 0;

  unsigned int diagramOrderHistogramLength = 1000;
  vector<unsigned long long int> diagramOrderHistogram(diagramOrderHistogramLength, 0);

  double secsSpentDoingDMC;
  do {



    for (unsigned int i = 0; i != numItrsPerTimesCheck; i++) {

      auto updateMethod = chooseUpdateMethod[this->Uint(0, chooseUpdateMethod.size() - 1)];
      (this->*updateMethod)(this->param);

      // bin diagrams of desired order and desired structure
      if (this->FD.Ds.size() >= this->minDiagramOrder) {

        if (this->FD.Ds.size() < diagramOrderHistogramLength) {
          diagramOrderHistogram[this->FD.Ds.size()]++;
        }

        meanTime += this->FD.length;
        totMeanTimes++;


        if (
          this->reducibleDiagrams ||
          ( this->skeletonDiagrams && this->FD.isSkeletonDiagram() ) ||
          ( ! this->skeletonDiagrams && this->FD.isIrreducibleDiagram() )
        ) {
          unsigned int ti = this->FD.length/this->dt;
          
          if (this->FD.Ds.size() != 1 || this->MCvsDMCboundary <= ti) {
            if (this->fixedExternalMomentum) {
              this->hist(0, ti)++;
            } else {
              unsigned int pi = this->FD.externalMomentum/this->dp;
              this->hist(pi, ti)++;
            }
          }
        }
      }

      // bin zeroth order diagrams used for normalization
      if (this->FD.Ds.size() == 0) {
        this->N0++;
      }
    }

    // increment current iteration number
    itrNum += numItrsPerTimesCheck;

    // time spent doing the DMC calculation
    secsSpentDoingDMC = double(clock() - beginTime)/CLOCKS_PER_SEC;

    // proposed new temporary save 
    unsigned int propTempSaveNum = secsSpentDoingDMC/secsPerTempSave;
    if (propTempSaveNum > tempSaveNum && propTempSaveNum < this->numTempDMCsaves + 1) {
      tempSaveNum = propTempSaveNum;

      if (true) {
        // time_t rawTime;
        // time(&rawTime);
        time_t rawTime = std::time(nullptr);

        char buffer[80];
        strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", localtime(&rawTime));
        string dateAndTimeString(buffer);


        cout << this->worldRank << " @ " << itrNum << " (" << secsSpentDoingDMC << ")" << endl;
      }

      cout << "[" << this->worldRank << "]" << endl;
      cout << "meanT_" << this->worldRank << " = " << meanTime/totMeanTimes << endl;
      cout << "\tNo_" << this->worldRank << " = np.array([";
      for (unsigned int i = 0; i < diagramOrderHistogramLength; i++) {
        cout << diagramOrderHistogram[i];
        if (i < diagramOrderHistogramLength - 1) {
          cout << ", ";
        }
      }
      cout << "])" << endl;
      cout << "\tN0_" << this->worldRank << " = " << this->N0 << endl;
      cout << "\tN_" << this->worldRank << " = np.array([";
      for (unsigned int i = 0; i < this->hist.size(); i++) {
        cout << this->hist(0, i);
        if (i < this->hist.size() - 1) {
          cout << ", ";
        }
      }
      cout << "])" << endl;

      if (this->worldSize > 1) {
        // sum up the contributions from each an every process
        Array<unsigned long long int, Dynamic, Dynamic> totHist;
        unsigned long long int totN0, totItrNum;
        this->sumHistograms(itrNum, totHist, totN0, totItrNum);

        if (this->worldRank == 0) {
          // write to file using totN0 and totHist
          this->write2file(totHist, totN0, totItrNum, secsSpentDoingDMC);
        }

      } else {
        // write to file using N0 and hist
        this->write2file(this->hist, this->N0, itrNum, secsSpentDoingDMC);
      }

      if (true) {
        // time_t rawTime;
        // time(&rawTime);
        time_t rawTime = std::time(nullptr);

        char buffer[80];
        strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", localtime(&rawTime));
        string dateAndTimeString(buffer);


        cout << this->worldRank << " @ cont " << " (" << secsSpentDoingDMC << ")" << endl;
      }

    }





    // // temporary saves
    // if (numItrAfterPrevTempSave++ == numItrBetweenTempSaves) {
    //   numItrAfterPrevTempSave = 0;

    //   if (true) {
    //     // time_t rawTime;
    //     // time(&rawTime);
    //     time_t rawTime = std::time(nullptr);

    //     char buffer[80];
    //     strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", localtime(&rawTime));
    //     string dateAndTimeString(buffer);


    //     cout << this->worldRank << " @ " << itrNum << " (" << dateAndTimeString << ")" << endl;
    //   }

    //   if (this->worldSize > 1) {
    //     // sum up the contributions from each an every process
    //     Array<unsigned long long int, Dynamic, Dynamic> totHist;
    //     unsigned long long int totN0, totItrNum;
    //     this->sumHistograms(totHist, totN0, totItrNum);

    //     if (this->worldRank == 0) {
    //       // write to file using totN0 and totHist
    //       this->write2file(totHist, totN0, totItrNum);
    //     }

    //   } else {
    //     // write to file using N0 and hist
    //     this->write2file(this->hist, N0, itrNum);
    //   }
    // }








  } while (secsSpentDoingDMC < this->numSecsDoingDMCperCorePerBoldItr);






  cout << "[Finished @ " << this->worldRank << " in "  << difftime(time(NULL), startTime) << "s]" << endl;

  if (this->worldSize > 1) {
    // sum up the contributions from each an every process
    Array<unsigned long long int, Dynamic, Dynamic> totHist;
    unsigned long long int totN0, totItrNum;

    this->sumHistograms(itrNum, totHist, totN0, totItrNum);

    // if (this->bold && this->boldIteration < this->numBoldIterations) {
    //   // store self energy for future use
    //   normalizeHistogram(totHist, totN0, this->S);
    // }

    if (this->worldRank == 0) {
      // write to file using totN0 and totHist
      this->write2file(totHist, totN0, totItrNum, secsSpentDoingDMC);
    }
  } else {
    // if (this->bold && this->boldIteration < this->numBoldIterations) {
    //   // store self energy for future use
    //   normalizeHistogram(this->hist, this->N0, this->S);
    // }

    // write to file using N0 and hist
    this->write2file(this->hist, N0, itrNum, secsSpentDoingDMC);
  }
}