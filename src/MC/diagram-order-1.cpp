#include "MonteCarlo.h"

void MonteCarlo::diagramOrder1 (double length, unsigned int index) {
  double value = 0;
  Vector3d P0 = this->externalMomentum;

  for (long unsigned int i = 0; i != this->numIterations; i++) {

    double std = 0, t1 = 0, t2 = 0, wInv_ts = 0;
    if (this->externalLegs) {
      // sample times
      t1 = this->Udouble(0, length),
      t2 = this->Udouble(t1, length),
      wInv_ts = length * (length - t1);

      std = pow(t2 - t1, -0.5);
    } else {
      std = pow(length, -0.5);
    }

    // sample momentum
    double
      q = abs(this->Ndouble(std)),
      theta = this->Udouble(0, M_PI),
      phi = this->Udouble(0, 2*M_PI),
      wInv_Q = 2*M_PI*M_PI * sqrt(0.5*M_PI*std*std) * exp(0.5*pow(q/std, 2.0));

    Vector3d Q{
      q*sin(theta)*cos(phi),
      q*sin(theta)*sin(phi),
      q*cos(theta)
    };

    // chemical potential factor

    if (this->externalLegs) {
      // integrand value
      double integrand = this->G(P0, 0, t1)
                       * this->G(P0 - Q, t1, t2) * this->D(Q, theta, t1, t2)
                       * this->G(P0, t2, length);

      value += integrand * wInv_ts * wInv_Q;
    } else {
      // integrand value
      double integrand = this->G(P0 - Q, 0, length) * this->D(Q, theta, 0, length);

      value += integrand * wInv_Q;
    }
  }

  this->values[index] += value/this->numIterations;
}