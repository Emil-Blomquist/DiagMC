#include "DiagrammaticMonteCarlo.h"

DiagrammaticMonteCarlo::DiagrammaticMonteCarlo (
  char **argv,
  string pathToConfigFile
) {
  this->argv = argv;

  MPI_Comm_size(MPI_COMM_WORLD, &this->worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &this->worldRank);

  // seed random generator
  unsigned long long int seed = chrono::system_clock::now().time_since_epoch().count() / this->worldSize * (this->worldRank + 1);
  // this->mt.seed(seed1);
  this->pcg.seed(seed);

  this->parseConfig(pathToConfigFile);

  // to reach a random start connfiguration
  this->untilStart = 10000000;

  // unique file name
  this->param = this->Udouble(0, 1);

  // store time at which the calculation began
  char buffer[80];
  time_t rawtime = time(nullptr);
  strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", localtime(&rawtime));
  this->dateAndTimeString = buffer;

  // bins of same size 
  this->maxLength = ceil(this->maxLength/dt)*this->dt;
  this->maxMomenta = ceil(this->maxMomenta/dp)*this->dp;

  // in order to create the histogram
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