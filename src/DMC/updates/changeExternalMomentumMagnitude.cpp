#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::changeExternalMomentumMagnitude (double param) {
  if (this->FD.length == 0) {
    // if time difference is zero momentum does not matter...
    cout << "DMC::changeExternalMomentumMagnitude: FD.length = 0 -> return" << endl;
    return;
  }

  double
    r = this->Udouble(0, 1),
    std = (this->FD.length < pow(10.0, -10.0) ? 100000 : 1/sqrt(this->FD.length)),
    p = sqrt(2)*std * boost::math::erf_inv(r*boost::math::erf(sqrt(0.5)*this->maxMomenta/std)),
    wInvp = sqrt(0.5*M_PI)*std*boost::math::erf(sqrt(0.5)*this->maxMomenta/std)*exp(0.5*pow(p/std, 2.0));

  Vector3d
    P{0, 0, p},
    dP = P - this->FD.ExternalMomentum,
    meanP = this->calculateMeanP(this->FD.start, this->FD.end);

  double
    oldp = this->FD.externalMomentum,
    oldwInvp = sqrt(0.5*M_PI)*std*boost::math::erf(sqrt(0.5)*this->maxMomenta/std)*exp(0.5*pow(oldp/std, 2.0));

  // contribution from boldification
  double boldContribution = 1;
  if (this->bold && this->boldIteration > 0) {
    boldContribution = 0;

    shared_ptr<Vertex> v = this->FD.start;
    while (v != this->FD.end) {
      boldContribution += this->additionalPhase(v->G[1]->momentum + dP, v->G[1]->end->position - v->position)
                        - this->additionalPhase(v->G[1]);

      v = v->G[1]->end;
    }

    boldContribution = exp(boldContribution);
  }

  double
    exponent = -(dP.dot(meanP) + 0.5*dP.squaredNorm())*this->FD.length,
    a = wInvp/oldwInvp * boldContribution * exp(exponent);

  double oldVal = 0;
  if (this->debug) {
    oldVal = this->evaluateDiagram();
  }

  // accept or reject update
  bool accepted = false;
  if (a > this->Udouble(0, 1)) {

    // set new momentum
    this->FD.setExternalMomentum(P, p, dP);

    accepted = true;
  }

  if (this->debug) {
    double val = this->evaluateDiagram();

    if (accepted) {
      this->checkAcceptanceRatio(boldContribution * exp(exponent)/(val/oldVal), "changeExternalMomentumMagnitude");
    }

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::changeExternalMomentumMagnitude " << endl
           << "a=" << a << endl
           << "accepted=" << accepted << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "p=" << p << endl
           << "P=" << P.transpose() << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "changeExternalMomentumMagnitude: " << accepted << " " << a << " " << boldContribution * exp(exponent)/(val/oldVal) << endl;
    }
  }
}