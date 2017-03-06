#include "../DiagrammaticMonteCarloV2.h"

void DiagrammaticMonteCarloV2::changeInternalPhononMomentumMagnitude (double param) {
  // requirement to lower: must be at least of order 1
  if (this->FD.Ds.size() == 0) {
    return;
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

  // is always accepted
  this->FD.setInternalPhononMomentum(d, Q);

  if ( ! isfinite(Q[0]) || ! isfinite(Q[1]) || ! isfinite(Q[2])) {
    cout
      << "-------------------------" << endl
      << "DMC::changeInternalPhononMomentumMagnitude: nan encountered" << endl
      << "Q=" << Q.transpose() << endl
      << "q=" << q << endl
      << "erf_inv=" << boost::math::erf_inv(param3) << endl
      << "param1=" << param1 << endl
      << "param2=" << param2 << endl
      << "param3=" << param3 << endl
      << "dt=" << d->end->position - d->start->position << endl
      << "d->theta=" << d->theta << endl
      << "P0=" << P0.transpose() << endl
      << "p0=" << p0 << endl
      << "-------------------------" << endl;
  }

  if (this->debug) {
    double val = this->FD();

    double a;
    if (val == 0) {
      a = 0;
    } else if (oldVal == 0) {
      a = 1;
    } else {
      // a = exp(log(val) - log(oldVal) - pow(param1, 2) * (pow(oldq - param2, 2) - pow(q - param2, 2)));
      a = val/oldVal * exp(-pow(param1, 2) * (pow(oldq - param2, 2) - pow(q - param2, 2)));
    }

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::changeInternalPhononMomentumMagnitude " << endl
           << "a=" << a << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "exp=" << exp(-pow(param1, 2) * (pow(oldq - param2, 2) - pow(q - param2, 2))) << endl
           << "param1=" << param1 << endl
           << "param2=" << param2 << endl
           << "Q=" << Q.transpose() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "changeInternalPhononMomentumMagnitude " << a << endl;
    }
  }
}




// -----------
// Overflow at DiagrammaticMonteCarloV2::calculateQ
// Q=nan nan nan
// P0= -nan -nan -nan
// Ep= 0 0 1
// Eo1= 1 0 0
// Eo2= 0 1 0
// q= nan
// theta= 0.981302
// phi= 5.98027
// -----------
// -----------
// Overflow at DiagrammaticMonteCarloV2::calculateQ
// Q=nan nan nan
// P0= nan nan nan
// Ep= 0 0 1
// Eo1= 1 0 0
// Eo2= 0 1 0
// q= nan
// theta= 3.14159
// phi= 0.633581
// -----------