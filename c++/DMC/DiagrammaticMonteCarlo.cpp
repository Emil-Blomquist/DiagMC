#include "DiagrammaticMonteCarlo.h"

DiagrammaticMonteCarlo::DiagrammaticMonteCarlo (
  Vector3d P,
  double length,
  double alpha,
  double mu
) : FD{P, length, alpha, mu} {
  // seed random generator
  unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
  this->mt.seed(seed1);

  // will print overflow exceptions etc.
  this->debug = false;
  // will print acceptance ratios etc. (debug must be true)
  this->loud = false;

  this->maxOrder = 3;
}


vector<double> DiagrammaticMonteCarlo::run (int numIterations) {
   // Display display(&this->FD);

  // vector of pointers to member function of Phonon
  vector<double (DiagrammaticMonteCarlo::*)(void)> updateMethods = {
    &DiagrammaticMonteCarlo::shiftVertexPosition,
    &DiagrammaticMonteCarlo::swapPhononConnections,
    &DiagrammaticMonteCarlo::changeInternalPhononMomentumDirection,
    &DiagrammaticMonteCarlo::changeInternalPhononMomentumMagnitude,
    &DiagrammaticMonteCarlo::raiseOrder,
    &DiagrammaticMonteCarlo::lowerOrder
  };

  // random start configuration
  for (int i = 0; i < 100; ++i) {
    auto updateMethod = updateMethods[this->Uint(0, updateMethods.size() - 1)];
    (this->*updateMethod)();
  }

  // bins for counting
  vector<double> bins(this->maxOrder + 1, 0);

  for (int i = 0; i < numIterations; ++i) {
    // choose update operation on random
    auto updateMethod = updateMethods[this->Uint(0, updateMethods.size() - 1)];

    // save configuration
    this->FD.save();

    // update diagram
    double a = (this->*updateMethod)();

    if (a < this->Udouble(0, 1)) {
      // rejected update
      this->FD.revert();
    }
    // display.render();

    bins[this->FD.Ds.size()]++;
  }

  for(vector<double>::size_type i = 0; i != bins.size(); i++) {
    bins[i] /= numIterations;
  }

  return bins;
}














double DiagrammaticMonteCarlo::Udouble (double from, double to) {
  uniform_real_distribution<double> distribution(from, to);
  return distribution(this->mt);
}

int DiagrammaticMonteCarlo::Uint (int fromIncluded, int toIncluded) {
  uniform_int_distribution<double> distribution(fromIncluded, toIncluded);
  return distribution(this->mt);
}

Vector3d DiagrammaticMonteCarlo::calculateP0 (shared_ptr<Phonon> d) {
  double dt = d->end->position - d->start->position;

  Vector3d P0(0, 0, 0);

  // loop through electrons between phonon propagator ends
  shared_ptr<Electron> g = d->start->G[1];
  while (g->start != d->end) {
    P0 += g->momentum*(g->end->position - g->start->position);

    g = g->end->G[1];
  }

  // add own momentum as well
  P0 = P0/dt + d->momentum;

  return P0;
}

Vector3d DiagrammaticMonteCarlo::calculateQ (Vector3d P0, double q, double theta, double phi) {
  Vector3d tempVector, Ep, Eo1, Eo2, Qp, Qo;
  
  Ep = P0.normalized();
  if (isnan(Ep[0])) {
    // if P0 ≈ (0,0,0) we use that theta is the angle against the z-axis
    Ep = Vector3d(0, 0, 1);
    Eo1 = Vector3d(1, 0, 0);
    Eo2 = Vector3d(0, 1, 0);
  } else {
    // find a vector orthogonal to Ep
    tempVector = Ep.cross(Ep + Vector3d(1, 0, 0));
    Eo1 = tempVector.normalized();
    // cout << "TEMP " << Eo1.transpose() << endl;
    if (isnan(Eo1[0])) {
      Eo1  = Ep.cross(Ep + Vector3d(0, 1, 0)).normalized();
    }

    // span rest of R^3
    Eo2 = Ep.cross(Eo1);
  }

  // calculate the parallel and the orthogonal component of Q
  Qp = q*cos(theta)*Ep;
  Qo = (Eo1*cos(phi) + Eo2*sin(phi)) * q*sin(theta);

  return Qp + Qo;
}