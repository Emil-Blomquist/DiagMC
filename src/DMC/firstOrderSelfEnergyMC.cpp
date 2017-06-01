#include "DiagrammaticMonteCarlo.h"

double DiagrammaticMonteCarlo::firstOrderSelfEnergyMC (double t, Vector3d P) {
  double value = 0;

  for (long unsigned int i = 0; i != this->numMCIterations; i++) {

    // sample momentum
    double
      std = 1/sqrt(t),
      q = abs(this->Ndouble(std)),
      theta = this->Udouble(0, M_PI),
      phi = this->Udouble(0, 2*M_PI),
      wInv_Q = 2*M_PI*M_PI * sqrt(0.5*M_PI*std*std) * exp(0.5*pow(q/std, 2.0));

    Vector3d Q{
      q*sin(theta)*cos(phi),
      q*sin(theta)*sin(phi),
      q*cos(theta)
    };

    double
      pq = (P - Q).norm(),
      omega = FeynmanDiagram::phononEnergy(q);
    
    // integrand value
    double integrand = Electron::value(pq, t, this->mu, this->dE(pq, t))
                     * Phonon::value(omega, t, theta, this->alpha);

    value += integrand * wInv_Q;
  }

  return value/this->numMCIterations;
}