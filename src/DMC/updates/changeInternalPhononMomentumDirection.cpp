#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::changeInternalPhononMomentumDirection (double param) {
  // requirement to lower: must be at least of order 1
  if (this->FD.Ds.size() == 0) {
    return;
  }

  double oldVal = 0;
  if (this->debug) {
    oldVal = this->evaluateDiagram();
  }

  // select internal phonon on random
  shared_ptr<Phonon> d = this->FD.Ds[this->Uint(0, this->FD.Ds.size() - 1)];

  double
    theta = this->Udouble(0, M_PI),
    oldTheta = d->theta,
    phi = this->Udouble(0, 2*M_PI),
    q = d->q;

  Vector3d
    Q{
      q*sin(theta)*cos(phi),
      q*sin(theta)*sin(phi),
      q*cos(theta)
    },
    oldQ = d->momentum,
    dQ = Q - oldQ,
    meanP = this->calculateMeanP(d->start, d->end);

  double
    sinOldTheta = sin(oldTheta),
    dt = d->end->position - d->start->position,
    dq2 = dQ.squaredNorm(),
    exponent = (dQ.dot(meanP) - 0.5*dq2)*dt;

  // contribution from boldification
  double boldContribution = 0;
  if (this->bold && this->boldIteration > 0) {
    shared_ptr<Vertex> v = d->start;
    while (v != d->end) {
      boldContribution += this->additionalPhase(v->G[1]->momentum - dQ, v->G[1]->end->position - v->position)
                        - this->additionalPhase(v->G[1]);
      
      v = v->G[1]->end;
    }
  }

  double a;
  if (sinOldTheta == 0) {
    a = 1;
  } else {
    a = sin(theta)/sinOldTheta * exp(exponent + boldContribution);
  }

  // accept or reject update
  bool accepted = false;
  if (a > this->Udouble(0, 1)) {

    // set new momentum
    this->FD.setInternalPhononMomentum(d, Q, q);

    // set angles corresponding to new momentum
    this->FD.setInternalPhononMomentumDirection(d, theta, phi);

    accepted = true;
  }

  if (this->debug) {
    double val = this->evaluateDiagram();

    if (accepted) {
      this->checkAcceptanceRatio(a / (val/oldVal), "changeInternalPhononMomentumDirection");
    }

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::changeInternalPhononMomentumDirection " << endl
           << "a=" << a << endl
           << "accepted=" << accepted << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "sin(oldTheta)=" << sin(oldTheta) << endl
           << "sin(theta)=" << sin(theta) << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "changeInternalPhononMomentumDirection: " << accepted << " " << a << " " << a / (val/oldVal) << endl;
    }
  }
}