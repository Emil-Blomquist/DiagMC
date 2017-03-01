#include "../DiagrammaticMonteCarloV1.h"

double DiagrammaticMonteCarloV1::changeInternalPhononMomentumMagnitude (double param) {
  // requirement to lower: must be at least of order 1
  if (this->FD.Ds.size() == 0) {
    return 0;
  }

  // select internal phonon on random
  shared_ptr<Phonon> d = this->FD.Ds[this->Uint(0, this->FD.Ds.size() - 1)];

  Vector3d P0 = this->calculateP0(d);
  double p0 = P0.norm();

  double
    param1 = sqrt(0.5*(d->end->position - d->start->position)),
    param2 = p0*cos(d->theta);

  double q,
    r = this->Udouble(0, 1),
    param3 = r + (r - 1)*erf(param1*param2);
  
  if (abs(param3) >= 1) {
    if (this->debug) {
      cout << "--------------------------------------------------------------------" << endl
           << "boost::math::erf_inv: Overflow Error prevented" << endl
           << "param1=" << param1 << endl
           << "param2=" << param2 << endl
           << "param3=" << param3 << endl
           << "--------------------------------------------------------------------" << endl;
    }

    param3 += (param3 > 0 ? -1 : 1)*DBL_EPSILON;
  }

  q = param2 + boost::math::erf_inv(param3)/param1;

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
    double
      val = this->FD(),
      a = exp(log(val) - log(oldVal) - pow(param1, 2) * (pow(oldq - param2, 2) - pow(q - param2, 2)));

    if (! ::isnan(a) && a < numeric_limits<double>::max()) {
      if (this->loud) { cout << "changeInternalPhononMomentumMagnitude " << a << endl; }
    } else {
      a = 0;
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::changeInternalPhononMomentumMagnitude" << endl
           << "a=" << a << endl
           << "param2=" << param2 << endl
           << "Q=" << Q.transpose() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "--------------------------------------------------------------------" << endl;
    }
  }
  
  return 1;
}