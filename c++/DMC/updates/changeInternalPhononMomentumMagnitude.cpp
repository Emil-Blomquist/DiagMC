#include "../DiagrammaticMonteCarlo.h"

double DiagrammaticMonteCarlo::changeInternalPhononMomentumMagnitude () {
  // select internal phonon on random
  shared_ptr<Phonon> d = this->FD.Ds[this->Uint(0, this->FD.Ds.size() - 1)];

  Vector3d P0 = this->calculateP0(d);
  double p0 = P0.norm();

  double
    param1 = sqrt(0.5*(d->end->position - d->start->position)),
    param2 = p0*cos(d->theta);

  double
    r = this->Udouble(0, 1),
    q = param2 + boost::math::erf_inv(r + (r - 1)*erf(param1*param2))/param2;

  Vector3d Q = this->calculateQ(P0, q, d->theta, d->phi);

  double oldq, oldVal;
  Vector3d oldQ;
  if (this->debug) {
    oldq = d->momentum.norm();
    oldQ = this->calculateQ(P0, oldq, d->theta, d->phi);
    this->FD.setInternalPhononMomentum(d, oldQ);
    oldVal = this->FD();
  }

  // update diagram
  this->FD.setInternalPhononMomentum(d, Q);

  if (this->debug) {
    double a, val = this->FD();

    a = val/oldVal;
    a *= exp(-pow(param1, 2) * (pow(oldq - param2, 2) - pow(q - param2, 2)));

    cout << "changeInternalPhononMomentumMagnitude: " << a << endl;
  }
  
  return 1;
}