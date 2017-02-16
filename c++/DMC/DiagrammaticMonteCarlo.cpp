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

  // will print acceptance ratios etc.
  this->debug = true;

  this->maxOrder = 3;

  //
  // TEMP
  //
  auto v1 = FD.insertVertex(0, 0.1);
  auto v2 = FD.insertVertex(1, 0.1);
  auto v3 = FD.insertVertex(2, 0.1);
  auto v4 = FD.insertVertex(3, 0.1);

  auto d1 = FD.addInternalPhonon(v1, v3, Vector3d(1,2,1), 1, 2);
  auto d2 = FD.addInternalPhonon(v2, v4, Vector3d(3,4,5), 1, 2);



  //
  // run the algorithm
  //

  Display display(&this->FD);


  // vector of pointers to member function of Phonon
  vector<double (DiagrammaticMonteCarlo::*)(void)> updateMethods = {
    // &DiagrammaticMonteCarlo::shiftVertexPosition,
    // &DiagrammaticMonteCarlo::swapPhononConnections,
    // &DiagrammaticMonteCarlo::changeInternalPhononMomentumDirection,
    &DiagrammaticMonteCarlo::changeInternalPhononMomentumMagnitude,
    // &DiagrammaticMonteCarlo::raiseOrder,
    // &DiagrammaticMonteCarlo::lowerOrder
  };


  // save configuration
  this->FD.save();

  double oldVal = this->FD();

  for (int i = 0; i < 100; ++i) {
    auto updateMethod = updateMethods[this->Uint(0, updateMethods.size() - 1)];
    // update diagram
    (this->*updateMethod)();
    display.render();
  }


  // revert to saved configuration
  this->FD.revert();

  double val = this->FD();

  cout << (val == oldVal) << endl;


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
  if (Ep == P0) {
    // if P0 ≈ (0,0,0) we use that theta is the angle against the z-axis
    Ep = Vector3d(0, 0, 1);
  }

  // find a vector orthogonal to Ep
  tempVector = Ep.cross(Ep + Vector3d(1, 0, 0));
  Eo1 = tempVector.normalized();
  if (tempVector == Eo1) {
    Eo1  = Ep.cross(Ep + Vector3d(0, 1, 0)).normalized();
  }

  // span rest of R^3
  Eo2 = Ep.cross(Eo1);

  // calculate the parallel and the orthogonal component of Q
  Qp = q*cos(theta)*Ep;
  Qo = (Eo1*cos(phi) + Eo2*sin(phi)) * q*sin(theta);

  return Qp + Qo;
}