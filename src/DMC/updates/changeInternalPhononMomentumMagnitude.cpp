#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::changeInternalPhononMomentumMagnitude (double param) {
  // requirement to lower: must be at least of order 1
  if (this->FD.Ds.size() == 0) {
    return;
  }

  // select internal phonon on random
  shared_ptr<Phonon> d = this->FD.Ds[this->Uint(0, this->FD.Ds.size() - 1)];

  if (d->end->position - d->start->position == 0) {
    // if time difference is zero momentum does not matter...
    cout << "DMC::changeInternalPhononMomentumMagnitude: dt = 0 -> return" << endl;
    return;
  }

  double oldVal = 0;
  if (this->debug) {
    oldVal = this->evaluateDiagram();
  }

  double
    t1 = d->start->position,
    t2 = d->end->position,
    std = (t2 - t1 < pow(10.0, -10.0) ? 100000 : 1/sqrt(t2 - t1)); 
  normal_distribution<double> normal(0.0, std);
  
  double
    q = abs(normal(this->pcg)),
    wInvq = sqrt(0.5*M_PI*pow(std, 2.0)) * exp(0.5*pow(q/std, 2.0));

  Vector3d
    Q{
      q*sin(d->theta)*cos(d->phi),
      q*sin(d->theta)*sin(d->phi),
      q*cos(d->theta)
    },
    oldQ = d->momentum,
    dQ = Q - oldQ,
    meanP = this->calculateMeanP(d->start, d->end);

  double
    oldq = d->q,
    oldwInvq = sqrt(0.5*M_PI*pow(std, 2.0)) * exp(0.5*pow(oldq/std, 2.0));

  double
    dq2 = dQ.squaredNorm(),
    exponent = (-this->FD.phononEnergy(q) + this->FD.phononEnergy(d->q) + dQ.dot(meanP) - 0.5*dq2)*(t2 - t1);

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

  double a = wInvq/oldwInvq * exp(exponent + boldContribution);

  // accept or reject update
  bool accepted = false;
  if (a > this->Udouble(0, 1)) {

    // set new momentum
    this->FD.setInternalPhononMomentum(d, Q, q);

    accepted = true;
  }

  if (this->debug) {
    double val = this->evaluateDiagram();

    if (accepted) {
      this->checkAcceptanceRatio(exp(exponent + boldContribution)/(val/oldVal), "changeInternalPhononMomentumMagnitude");
    }

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::changeInternalPhononMomentumMagnitude " << endl
           << "a=" << a << endl
           << "accepted=" << accepted << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "Q=" << Q.transpose() << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "changeInternalPhononMomentumMagnitude: " << accepted << " " << a << " " << exp(exponent + boldContribution)/(val/oldVal) << endl;
    }
  }
}