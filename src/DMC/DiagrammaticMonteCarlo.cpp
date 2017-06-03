#include "DiagrammaticMonteCarlo.h"

DiagrammaticMonteCarlo::DiagrammaticMonteCarlo (
  Vector3d P,
  double maxLength,
  double alpha,
  double mu,
  unsigned long int numIterations,
  // unsigned int numBins,
  double param,
  char **argv
) {
  // seed random generator
  unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
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
  this->fixedExternalMomentum = false;
  // if we want to use Dyson equation (fixedExternalMomentum must be false)
  this->Dyson = false;
  // wether or not we should employ boldification (fixedExternalMomentum must be false and Dyson must me set true)
  this->bold = true;

  // for when to bin the diagram
  this->minDiagramOrder = 1;
  // raise diagrm order will look at this (zero -> turned off)
  this->maxDiagramOrder = 5;
  // how many iterations in the bold scheme shall be done
  this->numBoldIterations = 0;

  // number of iterations used for each MC calculation
  this->numMCIterations = 10000000;

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

  // if we want to import a G this must be done before we initialize the Feynman diagram
  // this->importG("Glarge.txt");

  // initiate dE and dG (if needed) which are going to be overwritten for each iteration
  if ( ! this->fixedExternalMomentum) {
    const unsigned int
      Np = this->maxMomenta/this->dp,
      Nt = this->maxLength/this->dt;

    if (this->Dyson && ! this->dG.size()) {
      this->dG = Array<double, Dynamic, Dynamic>::Zero(Np, Nt);
    }
    
    if (this->Dyson && this->bold && ! this->dE.size()) {
      this->dE = Array<double, Dynamic, Dynamic>::Zero(Np, Nt);
    }
  }

  if (this->bold) {
    for (this->boldIteration = 0; this->boldIteration != this->numBoldIterations + 1; this->boldIteration++) {
      this->run();

      this->calculateEnergyDiff();


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
  const unsigned int Nt = this->maxLength/this->dt;
  if ( ! this->fixedExternalMomentum) {
    const unsigned int Np = this->maxMomenta/this->dp;

    this->hist = Array<unsigned long int, Dynamic, Dynamic>::Zero(Np, Nt);
  } else {

    this->bins = Array<unsigned long int, Dynamic, 1>::Zero(Nt, 1);
  }

  // save an empty file
  // this->write2file();

  // MC instead of DMC for S1 at small times
  if (this->minDiagramOrder <= 1 && ! this->externalLegs) {
    this->firstOrderSelfEnergyMC();
  }

  // // to reach start connfiguration
  // const unsigned int untilStart = 10000000;

  // // to save data under the process
  // const unsigned int saveAfter = 500*1000000;

  // multimap<unsigned int, void (DiagrammaticMonteCarlo::*)(double)> updateMethods;

  // // specify the relative probability of choosing a specific update function
  // if (this->bold && this->boldIteration > 0) {
  //   updateMethods = {
  //     {1, &DiagrammaticMonteCarlo::shiftVertexPosition},
  //     {1, &DiagrammaticMonteCarlo::swapPhononConnections},
  //     {1, &DiagrammaticMonteCarlo::changeInternalPhononMomentumDirection},
  //     {1, &DiagrammaticMonteCarlo::changeInternalPhononMomentumMagnitude},
  //     {1, &DiagrammaticMonteCarlo::BOLDraiseOrder}, // <- These two must have the same probability
  //     {1, &DiagrammaticMonteCarlo::BOLDlowerOrder}, // <-
  //     {5, &DiagrammaticMonteCarlo::BOLDchangeDiagramLength},
  //     {10, &DiagrammaticMonteCarlo::BOLDchangeDiagramLengthComplex},
  //     {(this->fixedExternalMomentum ? 0 : 1), &DiagrammaticMonteCarlo::changeExternalMomentumMagnitude}
  //   };
  // } else {
  //   updateMethods = {
  //     {1, &DiagrammaticMonteCarlo::shiftVertexPosition},
  //     {1, &DiagrammaticMonteCarlo::swapPhononConnections},
  //     {1, &DiagrammaticMonteCarlo::changeInternalPhononMomentumDirection},
  //     {1, &DiagrammaticMonteCarlo::changeInternalPhononMomentumMagnitude},
  //     {1, &DiagrammaticMonteCarlo::raiseOrder}, // <- These two must have the same probability
  //     {1, &DiagrammaticMonteCarlo::lowerOrder}, // <-
  //     {1, &DiagrammaticMonteCarlo::changeDiagramLength},
  //     {1, &DiagrammaticMonteCarlo::changeDiagramLengthComplex},
  //     {(this->fixedExternalMomentum ? 0 : 1), &DiagrammaticMonteCarlo::changeExternalMomentumMagnitude}
  //   };
  // }

  // // vector which is going to contain the specified quantity of update functions
  // vector<void (DiagrammaticMonteCarlo::*)(double)> chooseUpdateMethod;

  // // populate vector
  // for (auto updateMethod = updateMethods.begin(); updateMethod != updateMethods.end(); updateMethod++) {
  //   for (unsigned int i = 0; i != updateMethod->first; i++) {
  //     chooseUpdateMethod.push_back(updateMethod->second);
  //   }
  // }

  // // to start at a random position
  // for (unsigned int i = 0; i < untilStart; ++i) {
  //   // choose update operation on random
  //   auto updateMethod = chooseUpdateMethod[this->Uint(0, chooseUpdateMethod.size() - 1)];
  //   // update diagram
  //   (this->*updateMethod)(this->param);
  // }

  // // main loop
  // for (unsigned long int i = 0; i < this->numIterations; ++i) {
  //   auto updateMethod = chooseUpdateMethod[this->Uint(0, chooseUpdateMethod.size() - 1)];
  //   (this->*updateMethod)(this->param);

  //   if (i%saveAfter == saveAfter - 1) {
  //     // save temporary result
  //     this->write2file(i + 1);
  //   }

  //   // bin diagrams of desired order and desired structure
  //   if (this->FD.Ds.size() >= this->minDiagramOrder) {
  //     if (
  //       this->reducibleDiagrams ||
  //       ( this->skeletonDiagrams && this->FD.isSkeletonDiagram() ) ||
  //       ( ! this->skeletonDiagrams && this->FD.isIrreducibleDiagram() )
  //     ) {
  //       unsigned int ti = this->FD.length/this->dt;
        
  //       if (this->FD.Ds.size() != 1 || this->MCvsDMCboundary <= ti) {
  //         if (this->fixedExternalMomentum) {
  //           this->bins[ti]++;
  //         } else {
  //           unsigned int pi = this->FD.externalMomentum/this->dp;
  //           this->hist(pi, ti)++;
  //         }
  //       }
  //     }
  //   }

  //   // bin zeroth order diagrams used for normalization
  //   if (this->FD.Ds.size() == 0) {
  //     this->N0++;
  //   }
  // }

  // // save final result
  // this->write2file();
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