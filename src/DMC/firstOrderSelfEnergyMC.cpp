#include "DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::appendKeyValue (vector<KeyValue>& vec, unsigned int n) {
  unsigned int 
    pi = n/this->MCvsDMCboundary,
    ti = n - pi*this->MCvsDMCboundary;

  double
    p = (this->fixedExternalMomentum ? this->initialExternalMomentum.norm() : (pi + 0.5)*this->dp),
    t = (ti + 0.5)*this->dt,
    S1 = firstOrderSelfEnergyMC(t, Vector3d{0, 0, p});

  KeyValue keyValue {pi, ti, S1};
  vec.push_back(keyValue);
}

void DiagrammaticMonteCarlo::firstOrderSelfEnergyMC () {

  // to time the time needed
  clock_t tStart = clock();

  // calculate MC/DMC boundary for S1
  // doing this for p = dG = 0 so we get an indication
  double threshold = pow(10, -3.0);
  unsigned int i = 0;
  do {
    double
      t_l = i*this->dt,
      t = (i + 0.5)*this->dt,
      t_u = (i + 1)*this->dt,
      param1 = sqrt(FeynmanDiagram::phononEnergy(1) - this->mu),
      meanValue = this->alpha/(param1*this->dt) * (erf(param1*sqrt(t_u)) - erf(param1*sqrt(t_l))),
      trueValue = this->alpha/sqrt(M_PI*t) * exp(-pow(param1, 2.0)*t);

      if (meanValue - trueValue < threshold) {
        this->MCvsDMCboundary = i;
        break;
      }

  } while (++i);

  if (this->worldSize > 1) {
    // parallelize over all processes

    unsigned int
      Nt = this->MCvsDMCboundary,
      numEach = this->Np*Nt/this->worldSize,
      from = this->worldRank*numEach,
      end = (this->worldRank + 1)*numEach;

    vector<KeyValue> vector2send, vector2receive{this->Np*Nt};

    for (unsigned int n = from; n != end; n++) {
      this->appendKeyValue(vector2send, n);
    }

    // rest
    unsigned int myRest = this->worldSize*numEach + this->worldRank;
    if (myRest <= this->Np*Nt - 1) {
      this->appendKeyValue(vector2send, myRest);
    }

    // calculate receive counts and corresponding displacements
    vector<int> recvCounts(this->worldSize), recvDiscps(this->worldSize);
    if (this->worldRank == 0) {
      for (int rank = 0; rank != this->worldSize; rank++) {
        int count = (numEach + (this->worldSize*numEach + rank <= this->Np*Nt - 1 ? 1 : 0)) * sizeof(KeyValue);

        recvCounts[rank] = count;
        recvDiscps[rank] = (rank == 0 ? 0 : recvDiscps[rank - 1] + recvCounts[rank - 1]);
      }
    }

    // Gather all data at root
    MPI_Gatherv(
      vector2send.data(),
      vector2send.size() * sizeof(KeyValue),
      MPI_BYTE,
      vector2receive.data(),
      &recvCounts[0],
      &recvDiscps[0],
      MPI_BYTE,
      0,
      MPI_COMM_WORLD);

    // store the S1
    if (this->worldRank == 0) {
      this->S1mc = ArrayXXd::Zero(this->Np, Nt);
      for (KeyValue kv : vector2receive) {
        this->S1mc(kv.pi, kv.ti) = kv.val;
      }
    }

    printf("[S1 MC @ %u in %.2fs]\n", this->worldRank, (double)(clock() - tStart)/CLOCKS_PER_SEC);
  } else {
    // MC calculation
    this->S1mc = Array<double, Dynamic, Dynamic>::Zero(this->Np, this->MCvsDMCboundary);

    for (unsigned int i = 0; i != this->S1mc.rows(); i++) {
      double p = (this->fixedExternalMomentum ? this->initialExternalMomentum.norm() : (i + 0.5)*this->dp);
      for (unsigned int j = 0; j != this->S1mc.cols(); j++) {
        double t = (j + 0.5)*this->dt;

        this->S1mc(i, j) = firstOrderSelfEnergyMC(t, Vector3d{0, 0, p});
      }
    }

    printf("[S1 MC in %.2fs]\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  }
}

double DiagrammaticMonteCarlo::firstOrderSelfEnergyMC (double t, Vector3d P) {
  double value = 0;

  for (long unsigned int i = 0; i != this->numMCIterations; i++) {
    // sample momentum
    double
      std = 1/sqrt(t),
      q = abs(this->Ndouble(std)),
      theta = this->Udouble(0, M_PI),
      phi = this->Udouble(0, 2*M_PI),
      wInv_Q = 2*M_PI*M_PI * sqrt(0.5*M_PI*std*std) * exp(0.5*pow(q/std, 2.0));

    Vector3d Q{
      q*sin(theta)*cos(phi),
      q*sin(theta)*sin(phi),
      q*cos(theta)
    };

    double
      pq = (P - Q).norm(),
      omega = FeynmanDiagram::phononEnergy(q);
    
    // integrand value
    double integrand = Electron::value(pq, t, this->mu, (this->dE.size() ? this->dEOf(pq, t) : 0))
                     * Phonon::value(omega, t, theta, this->alpha);

    value += integrand * wInv_Q;
  }

  return value/this->numMCIterations;
}