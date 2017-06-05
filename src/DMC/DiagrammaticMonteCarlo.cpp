#include "DiagrammaticMonteCarlo.h"

DiagrammaticMonteCarlo::DiagrammaticMonteCarlo (
  Vector3d P,
  double maxLength,
  double alpha,
  double mu,
  unsigned long long int numIterations,
  // unsigned int numBins,
  double param,
  char **argv
) {
  MPI_Comm_size(MPI_COMM_WORLD, &this->worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &this->worldRank);

  // seed random generator
  unsigned long long int seed1 = chrono::system_clock::now().time_since_epoch().count() / this->worldSize * (this->worldRank + 1);
  this->mt.seed(seed1);

  // will print overflow exceptions etc. (code is half as fast with this activated)
  this->debug = false;
  // will print acceptance ratios etc. (debug must be true)
  this->loud = false;

  // if the diagram should have external legs or not
  this->externalLegs = false;
  // if we should only take into account irreducible diagrams
  this->reducibleDiagrams = false;
  // if we should sample only skeleton diagrams (reducibleDiagrams must be false)
  this->skeletonDiagrams = true;
  // if we should let the external momentum vary or not
  this->fixedExternalMomentum = true;
  // if we want to use Dyson equation (fixedExternalMomentum must be false)
  this->Dyson = false;
  // wether or not we should employ boldification (fixedExternalMomentum must be false and Dyson must me set true)
  this->bold = true;

  // for when to bin the diagram
  this->minDiagramOrder = 1;
  // raise diagrm order will look at this (zero -> turned off)
  this->maxDiagramOrder = 10;
  // how many iterations in the bold scheme shall be done
  this->numBoldIterations = 1;

  // number of iterations used for each MC calculation
  this->numMCIterations = 1000000;
  // to reach a random start connfiguration
  this->untilStart = 10000000;
  // how often in seconds we should save by writing to file
  this->savePeriod = 10;

  this->initialExternalMomentum = P;
  this->mu = mu;
  this->alpha = alpha;
  this->numIterations = numIterations;
  this->param = this->Udouble(0, 1);
  this->argv = argv;

  // store time at which the calculation began
  time_t rawtime;
  time (&rawtime);
  this->timeinfo = localtime(&rawtime);

  // in order to create the histogram
  this->maxLength = maxLength;
  this->dt = 0.02;

  this->maxMomenta = 0.08;
  this->dp = 0.02;

  this->Np = (this->fixedExternalMomentum ? 1 : this->maxMomenta/this->dp);
  this->Nt = this->maxLength/this->dt;

  // if we want to import a G this must be done before we initialize the Feynman diagram
  this->importG("Glarge.txt");

  if (this->bold) {
    for (this->boldIteration = (this->dG.size() ? 1 : 0); this->boldIteration != this->numBoldIterations + 1; this->boldIteration++) {
      this->run();

      if (this->boldIteration < this->numBoldIterations) {
        // perform Dyson using this->S
        this->doDyson(this->S, this->dG);

        // calculate the energy difference and store in this->dE
        this->calculateEnergyDiff(dG, this->dE);
      }



      // Array<double, Dynamic, Dynamic> normHist = ArrayXXd::Zero(this->hist.rows(), this->hist.cols());
      // this->normalizedHistogram(normHist);


      // // normalized histogram 
      // cout << "normHist" << this->boldIteration + 1 << " = np.array([";
      // for (unsigned int i = 0; i != normHist.cols(); i++) {
      //   cout << normHist(0, i);

      //   if (i < normHist.cols() - 1) {
      //     cout << ", ";
      //   }
      // }
      // cout << "])" << endl << endl;


      // cout << "dG" << this->boldIteration + 1 << " = np.array([";
      // for (unsigned int i = 0; i != this->dG.cols(); i++) {
      //   cout << this->dG(0, i);

      //   if (i < this->dG.cols() - 1) {
      //     cout << ", ";
      //   }
      // }
      // cout << "])" << endl << endl;

      // cout << "---------" << endl;

      // bold iteration is complete, compute dG to be used for future histogram normalizations
      // this->doDyson(this->dG);
    }

  } else {
    this->run();
  }
}

void DiagrammaticMonteCarlo::run () {
  // initiate/reset FeynmanDiagram
  this->FD = FeynmanDiagram{this->initialExternalMomentum, this->dt, this->alpha, this->mu};

  // initiate/reset histogram
  this->N0 = 0;
  this->hist = Array<unsigned long long int, Dynamic, Dynamic>::Zero(this->Np, this->Nt);

  // save an empty file
  if (this->worldRank == 0) {
    this->write2file(this->hist, this->N0, this->currItr);
  }

  // used for temporary saving result
  time_t startTime = time(NULL);

  // MC instead of DMC for S1 at small times
  if (this->minDiagramOrder <= 1 && ! this->externalLegs) {
    this->firstOrderSelfEnergyMC();
  }

  // specify the relative probability of choosing a specific update function
  multimap<unsigned int, void (DiagrammaticMonteCarlo::*)(double)> updateMethods;
  if (this->bold && this->boldIteration > 0) {
    updateMethods = {
      {1, &DiagrammaticMonteCarlo::shiftVertexPosition},
      {1, &DiagrammaticMonteCarlo::swapPhononConnections},
      {1, &DiagrammaticMonteCarlo::changeInternalPhononMomentumDirection},
      {1, &DiagrammaticMonteCarlo::changeInternalPhononMomentumMagnitude},
      {1, &DiagrammaticMonteCarlo::BOLDraiseOrder}, // <- These two must have the same probability
      {1, &DiagrammaticMonteCarlo::BOLDlowerOrder}, // <-
      {5, &DiagrammaticMonteCarlo::BOLDchangeDiagramLength},
      {10, &DiagrammaticMonteCarlo::BOLDchangeDiagramLengthComplex},
      {(this->fixedExternalMomentum ? 0 : 1), &DiagrammaticMonteCarlo::changeExternalMomentumMagnitude}
    };
  } else {
    updateMethods = {
      {1, &DiagrammaticMonteCarlo::shiftVertexPosition},
      {1, &DiagrammaticMonteCarlo::swapPhononConnections},
      {1, &DiagrammaticMonteCarlo::changeInternalPhononMomentumDirection},
      {1, &DiagrammaticMonteCarlo::changeInternalPhononMomentumMagnitude},
      {1, &DiagrammaticMonteCarlo::raiseOrder}, // <- These two must have the same probability
      {1, &DiagrammaticMonteCarlo::lowerOrder}, // <-
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
  
  // to reach a random start connfiguration
  for (unsigned long long int i = 0; i < this->untilStart; ++i) {
    // choose update operation on random
    auto updateMethod = chooseUpdateMethod[this->Uint(0, chooseUpdateMethod.size() - 1)];
    // update diagram
    (this->*updateMethod)(this->param);
  }

  // number of iterations for each process
  unsigned long long int
    localNumIterations = this->numIterations/this->worldSize,
    rest = this->numIterations - localNumIterations*this->worldSize;
  if (this->worldRank < rest) localNumIterations++;

  // main loop
  for (this->currItr = 0; this->currItr < localNumIterations; this->currItr++) {
    auto updateMethod = chooseUpdateMethod[this->Uint(0, chooseUpdateMethod.size() - 1)];
    (this->*updateMethod)(this->param);

    // bin diagrams of desired order and desired structure
    if (this->FD.Ds.size() >= this->minDiagramOrder) {
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

    // temporary saves
    if (difftime(time(NULL), startTime) > this->savePeriod) {
      if (this->worldSize > 1) {
        // sum up the contributions from each an every process
        Array<unsigned long long int, Dynamic, Dynamic> totHist;
        unsigned long long int totN0, totCurrItr;
        this->sumHistograms(totHist, totN0, totCurrItr);

        // when all processes are synchronized, reset the timer
        startTime = time(NULL);

        if (this->worldRank == 0) {
          // write to file using totN0 and totHist
          this->write2file(totHist, totN0, totCurrItr);
        }

      } else {
        // write to file using N0 and hist
        this->write2file(this->hist, N0, this->currItr);
        startTime = time(NULL);
      }
    }
  }

  if (this->worldSize > 1) {
    // sum up the contributions from each an every process
    Array<unsigned long long int, Dynamic, Dynamic> totHist;
    unsigned long long int totN0, totCurrItr;
    this->sumHistograms(totHist, totN0, totCurrItr);

    if (this->bold && this->boldIteration < this->numBoldIterations) {
      // store self energy for future use
      normalizeHistogram(totHist, totN0, this->S);
    }

    if (this->worldRank == 0) {
      // write to file using totN0 and totHist
      this->write2file(totHist, totN0, totCurrItr);
    }
  } else {
    if (this->bold && this->boldIteration < this->numBoldIterations) {
      // store self energy for future use
      normalizeHistogram(this->hist, this->N0, this->S);
    }

    // write to file using N0 and hist
    this->write2file(this->hist, N0, this->currItr);
  }
}








double DiagrammaticMonteCarlo::Udouble (double from, double to) {
  uniform_real_distribution<double> distribution(from, to);
  return distribution(this->mt);
}

double DiagrammaticMonteCarlo::Ndouble (double std) {
  normal_distribution<double> distribution(0.0, std);
  return distribution(this->mt);
}

int DiagrammaticMonteCarlo::Uint (int fromIncluded, int toIncluded) {
  uniform_int_distribution<int> distribution(fromIncluded, toIncluded);
  return distribution(this->mt);
}

Vector3d DiagrammaticMonteCarlo::calculateMeanP (shared_ptr<Vertex> v1, shared_ptr<Vertex> v2) {

  Vector3d meanP(0, 0, 0);

  if (v1 != v2) {
    if (v1->G[1]->end == v2) {
      // if dt = 0
      meanP += v1->G[1]->momentum;

      if ( ! isfinite(meanP[0]) || ! isfinite(meanP[1]) || ! isfinite(meanP[2])) {
        cout << "-----------" << endl
             << "Overflow at DiagrammaticMonteCarlo::calculateMeanP" << endl
             << "meanP= " << meanP.transpose() << endl
             << "t1= " << v1->position << endl
             << "t2= " << v2->position << endl
             << "-----------" << endl;
      }
    } else {
      double dt = v2->position - v1->position;

      // loop through electrons between phonon propagator ends
      shared_ptr<Electron> g = v1->G[1];
      do {
        meanP += g->momentum*(g->end->position - g->start->position);
      } while (g->end != v2 && (g = g->end->G[1]));

      if (dt == 0 || ! isfinite(meanP[0]) || ! isfinite(meanP[1]) || ! isfinite(meanP[2])) {
        cout << "-----------" << endl
             << "Overflow at DiagrammaticMonteCarlo::calculateMeanP" << endl
             << "dt=" << dt << endl
             << "meanP= " << meanP.transpose() << endl
             << "t1= " << v1->position << endl
             << "t2= " << v2->position << endl
             << "-----------" << endl;
      }

      meanP /= dt;
    }
  }

  return meanP;
}

Vector3d DiagrammaticMonteCarlo::calculateP0 (shared_ptr<Phonon> d) {
  Vector3d meanP = this->calculateMeanP(d->start, d->end);

  return meanP + d->momentum;
}

void DiagrammaticMonteCarlo::checkAcceptanceRatio (double a, string updateFunction) {
  if (abs(a - 1) > pow(10.0, -5)) {
    cout << "Acceptance ratio at " << updateFunction << ": "  << "n=" << this->FD.Ds.size() << " -> " << setprecision(17) << a - 1 << endl;
  }
}