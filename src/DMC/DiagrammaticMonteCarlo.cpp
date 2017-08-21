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
  unsigned long long int seed = chrono::system_clock::now().time_since_epoch().count() / this->worldSize * (this->worldRank + 1);
  // this->mt.seed(seed1);
  this->pcg.seed(seed);

  // will print overflow exceptions etc. (code is half as fast with this activated)
  this->debug = false;
  // will print acceptance ratios etc. (debug must be true)
  this->loud = false;

  // if the diagram should have external legs or not
  this->externalLegs = true;
  // if we should only take into account irreducible diagrams
  this->reducibleDiagrams = true;
  // if we should sample only skeleton diagrams (reducibleDiagrams must be false)
  this->skeletonDiagrams = false;
  // if we should let the external momentum vary or not
  this->fixedExternalMomentum = true;
  // if we want to use Dyson equation (fixedExternalMomentum must be false (REALLY!?))
  this->Dyson = false; // if we want to output using dyson or not
  // wether or not we should employ boldification (fixedExternalMomentum must be false and Dyson must me set true)
  this->bold = false;

  // for when to bin the diagram
  this->minDiagramOrder = 0;
  // raise diagrm order will look at this (zero -> turned off)
  this->maxDiagramOrder = 0;
  // how many iterations in the bold scheme shall be done
  this->numBoldIterations = 4;

  // number of iterations used for each MC calculation
  // this->numMCIterations = 100000;
  // to reach a random start connfiguration
  this->untilStart = 10000000;
  // how often in seconds we should save by writing to file
  this->numTempDMCsaves = 9;





  // numSecsPerCorePerBoldItr = 10 * 60;
  // fracToSpendOnMC = 0.1;

  this->numSecsDoingMCperCorePerBoldItr = 30; // 30 min
  this->numSecsDoingDMCperCorePerBoldItr = 10*60;//10*60*60; // 10 h





  this->initialExternalMomentum = P;
  this->mu = mu;
  this->alpha = alpha;
  // this->numIterations = numIterations;
  this->param = this->Udouble(0, 1);
  this->argv = argv;

  // store time at which the calculation began
  char buffer[80];
  time_t rawtime = time(nullptr);
  strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", localtime(&rawtime));
  this->dateAndTimeString = buffer;


  // in order to create the histogram
  this->maxLength = maxLength;
  this->dt = 0.01;

  this->maxMomenta = 1;
  this->dp = 0.02;

  this->Np = (this->fixedExternalMomentum ? 1 : this->maxMomenta/this->dp);
  this->Nt = this->maxLength/this->dt;

  // if we want to import a G this must be done before we initialize the Feynman diagram
  // this->importG("Glarge.txt");
  // this->importS("Slarge.txt");

  if (this->bold) {
    for (
      this->boldIteration = (this->dG.size() ? this->boldIteration + 1 : 0);
      this->boldIteration != this->numBoldIterations + 1;
      this->boldIteration++
    ) {
      this->run();

      if (this->boldIteration < this->numBoldIterations) {

        if (this->worldRank == 0) {
          ArrayXXd S;
          this->normalizeHistogram(this->hist, this->N0, S);

          // perform Dyson using this->S
          this->doDyson(S, this->dG);
        } else {
          // for the other ranks, create a buffer large enough to contain the data
          this->dG = ArrayXXd::Zero(this->hist.rows(), this->hist.cols());
        }

        // broadcast dG to other ranks as well
        MPI_Bcast(
          this->dG.data(),
          this->dG.size(),
          MPI_DOUBLE,
          0,
          MPI_Comm MPI_COMM_WORLD);

        // calculate the energy difference and store in this->dE
        this->calculateEnergyDiff(this->dG, this->dE);

        // calculate the rate parameters of the exponential imaginary-time distributions
        this->calculateLambdas(this->dE, this->lambdas);
      }
    }
  } else {
    this->run();
  }
}










double DiagrammaticMonteCarlo::Udouble (double from, double to) {
  uniform_real_distribution<double> distribution(from, to);
  // return distribution(this->mt);
  return distribution(this->pcg);
}

double DiagrammaticMonteCarlo::Ndouble (double std) {
  normal_distribution<double> distribution(0.0, std);
  // return distribution(this->mt);
  return distribution(this->pcg);
}

int DiagrammaticMonteCarlo::Uint (int fromIncluded, int toIncluded) {
  uniform_int_distribution<int> distribution(fromIncluded, toIncluded);
  // return distribution(this->mt);
  return distribution(this->pcg);
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