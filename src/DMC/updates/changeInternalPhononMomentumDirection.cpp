#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::changeInternalPhononMomentumDirection (double param) {
  // requirement to lower: must be at least of order 1
  if (this->FD.Ds.size() == 0) {
    return;
  }

  // select internal phonon on random
  shared_ptr<Phonon> d = this->FD.Ds[this->Uint(0, this->FD.Ds.size() - 1)];

  double
    theta = this->Udouble(0, M_PI),
    oldTheta = d->theta,
    phi = this->Udouble(0, 2*M_PI),
    q = d->momentum.norm();

  Vector3d
    Q(q*sin(theta)*cos(phi), q*sin(theta)*sin(phi), q*cos(theta)),
    oldQ = d->momentum,
    dQ = Q - oldQ,
    meanP = this->calculateMeanP(d->start, d->end);


  double
    sinOldTheta = sin(oldTheta),
    dt = d->end->position - d->start->position,
    dq2 = dQ.squaredNorm(),
    exponent = (dQ.dot(meanP) - 0.5*dq2)*dt;


  double a;
  if (sinOldTheta == 0) {
    a = 1;
  } else {
    a = sin(theta)/sinOldTheta * exp(exponent);
  }


  double oldVal = 0;
  if (this->debug) {
    oldVal = this->FD();
  }

  // accept or reject update
  bool accepted = false;
  if (a > this->Udouble(0, 1)) {

    // set new momentum
    this->FD.setInternalPhononMomentum(d, Q);

    // set angles corresponding to new momentum
    this->FD.setInternalPhononMomentumDirection(d, theta, phi);

    accepted = true;
  }

  if (this->debug) {
    double
      val = this->FD(),
      acc = val/oldVal;

    // cout << "changeInternalPhononMomentumDirection " << (abs(a - acc) < pow(10.0, -10)) << endl;

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::changeInternalPhononMomentumDirection " << endl
           << "a=" << a << endl
           << "acc=" << acc << endl
           << "accepted=" << accepted << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "sin(oldTheta)=" << sin(oldTheta) << endl
           << "sin(theta)=" << sin(theta) << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "changeInternalPhononMomentumDirection " << a  << " " << acc << endl;
    }
  }
}