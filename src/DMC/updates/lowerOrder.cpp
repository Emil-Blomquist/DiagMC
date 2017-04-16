#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::lowerOrder (double param) {
  // requirement to lower: must be at least of order 1
  if (this->FD.Ds.size() == 0) {
    return;
  }

  // old configuration value
  double oldVal = 0;
  if (this->debug) {
    oldVal = this->FD();
  }

  // pick a phonon line
  int phononIndex = this->Uint(0, this->FD.Ds.size() - 1);
  shared_ptr<Phonon> d = this->FD.Ds[phononIndex];
  double Winvd = this->FD.Ds.size();

  // vertices
  shared_ptr<Vertex>
    v1 = d->start,
    v2 = d->end;

  double wInvQ, wInvt1, wInvt2, wInvG1, wInvd;

  if ( ! this->externalLegs && this->FD.Ds.size() == 1) {
    // go down to zeroth order

    wInvt1 = 1;
    wInvt2 = 1;
    wInvG1 = 1;
    wInvd = 1;

    double std = (v2->position - v1->position < pow(10.0, -10.0) ? 100000 : 1/sqrt(v2->position - v1->position)); 
    normal_distribution<double> normal(0.0, std);

    wInvQ = 2*pow(M_PI, 2.0) * sqrt(0.5*M_PI*pow(std, 2.0)) * exp(0.5*pow(d->q/std, 2.0));

  } else {
    if ( ! this->externalLegs && (v1 == this->FD.start || v2 == this->FD.end)) {
      // invalid diagram to lower
      return;
    }

    // calculate vertex probability
    wInvG1 = this->FD.Gs.size() - 2;

    // inverse probabilities for internal parameters
    // must be the same probabilities as for raising order!
    wInvt1 = (v1->G[1]->end == v2 ? v2->G[1]->end->position : v1->G[1]->end->position) - v1->G[0]->start->position;

    double
      t2low = v1->position,
      t2up = this->FD.length,
      l = 0.01,
      dt2 = v2->position - t2low;

    wInvt2 = exp(l*dt2)*(1 - exp(-l*(t2up - t2low)))/l;

    double std = (v2->position - v1->position < pow(10.0, -10.0) ? 100000 : 1/sqrt(v2->position - v1->position));

    wInvQ = 2*pow(M_PI, 2.0) * sqrt(0.5*M_PI*pow(std, 2.0)) * exp(0.5*pow(d->q/std, 2.0));
  }

  // to calculate the acceptance ratio
  Vector3d
    Q = d->momentum,
    P0 = this->calculateP0(d);

  double
    sinTheta = sin(d->theta),
    alpha = this->alpha,
    q2 = d->q*d->q,
    dt = v2->position - v1->position,
    exponential = exp(dt*(this->FD.phononEnergy(d->q) + 0.5*q2 - Q.dot(P0)));

  double a;
  if (wInvQ == 0 || sinTheta == 0 || wInvt2 == 0 || ! isfinite(exponential)) {
    a = 1;
  } else {
    a = exponential * sqrt(8)*M_PI*M_PI/(alpha*sinTheta) * Winvd/(wInvG1*wInvt1*wInvt2*wInvQ);
  }

  // accept or reject update
  bool accepted = false;
  if (a > this->Udouble(0, 1)) {
    // remove phonon
    this->FD.removeInternalPhonon(phononIndex);

    if (
      this->externalLegs ||
      ( ! this->externalLegs && this->FD.Ds.size() > 0))
    {
      // remove vertices
      this->FD.removeVertex(v1);
      this->FD.removeVertex(v2);
    }

    accepted = true;
  }

  if (this->debug) {
    double val = this->FD();

    if (accepted) {
      this->checkAcceptanceRatio(exponential * sqrt(8)*M_PI*M_PI/(alpha*sinTheta) /(val/oldVal), "lowerOrder");
    }

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::lowerOrder " << endl
           << "accepted=" << accepted << endl
           << "a=" << a << endl
           << "a_diag=" << val/oldVal * Winvd/(wInvG1*wInvt1*wInvt2*wInvQ) << endl
           << "ratio=" << a/(val/oldVal * Winvd/(wInvG1*wInvt1*wInvt2*wInvQ)) << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "wInvQ=" << wInvQ << endl
           << "wInvt2=" << wInvt2 << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "lowerOrder: " << accepted << " " << a << " " << exponential * sqrt(8)*M_PI*M_PI/(alpha*sinTheta) /(val/oldVal) << endl;
    }
  }
}