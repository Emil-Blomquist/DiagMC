#include "../DiagrammaticMonteCarlo.h"

double DiagrammaticMonteCarlo::changeInternalPhononMomentumDirection () {
  // requirement to lower: must be at least of order 1
  if (this->FD.Ds.size() == 0) {
    return 0;
  }

  // select internal phonon on random
  shared_ptr<Phonon> d = this->FD.Ds[this->Uint(0, this->FD.Ds.size() - 1)];

  Vector3d P0 = this->calculateP0(d);
  double
    q = d->momentum.norm(),
    p0 = P0.norm();

  double
    r = this->Udouble(0, 1),
    phi = this->Udouble(0, 2*M_PI);

  // sample new theta
  double cosTheta, param = q*p0*(d->end->position - d->start->position);
  if (param < pow(10.0, -10.0)) {
    // small "param" approximation to avoid overflow
    cosTheta = 1 - 2*r;
  } else {
    cosTheta = 1 + log(1 - r*(1 - exp(-2*param)))/param;
  }
  double theta = acos(cosTheta);

  // calculate new Q
  Vector3d Q = calculateQ(P0, q, theta, phi);

  double oldCosTheta, oldTheta, oldVal;
  if (this->debug) {
    // if P0 ≈ (0,0,0) we use that theta is the angle against the z-axis
    if (isnan(P0.normalized()[0])) {
      oldCosTheta = d->momentum[2]/q;
    } else {
      oldCosTheta = d->momentum.dot(P0)/(q*p0);
    }

    // deal with rounding errors
    if (oldCosTheta <= -1) {
      oldTheta = M_PI;
    } else {
      oldTheta = acos(oldCosTheta);
    }

    //
    // "d->theta" is outdated since the diagram might have changed in the sense that
    // it is no longer valid for comparing Q with a new Q'. Even the same Q might be
    // more/less preffered. "d->theta" is only good for evaluating the diagram (not this evaluation however)
    // -> we need to update the theta angle befor comparing with the new one (phi has not changed)
    //
    this->FD.setInternalPhononMomentumDirection(d, oldTheta, d->phi);

    // calculate the value correcponding to the old configuration
    oldVal = this->FD();
  }

  // update diagram (the last one only saves the angles)
  this->FD.setInternalPhononMomentum(d, Q);
  this->FD.setInternalPhononMomentumDirection(d, theta, phi);

  if (this->debug) {
    double val = this->FD();

    double a;
    if (abs(oldVal) > 0 and sin(theta) > 0) {
      a = val/oldVal;
      a *= sin(oldTheta)/sin(theta);
      a *= exp(-param*(cos(theta) - cos(oldTheta)));
      if (this->loud) { cout << "changeInternalPhononMomentumDirection " << a << endl; }
    } else {
      a = 0;
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::changeInternalPhononMomentumDirection Q=" << Q.transpose() << " oldVal=" << oldVal << " sinTheta=" << sin(theta) << endl
           << "--------------------------------------------------------------------" << endl;
    }
  }

  return 1;
}