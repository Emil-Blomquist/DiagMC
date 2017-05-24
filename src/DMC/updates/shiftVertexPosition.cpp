#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::shiftVertexPosition (double param) {
  // requirement to lower: must be at least of order 1
  if (
    this->FD.Ds.size() == 0 ||
    ( ! this->externalLegs && this->FD.Ds.size() == 1)
  ) {
    return;
  }

  // select vertex on random
  shared_ptr<Electron> g = this->FD.Gs[this->Uint(0, this->FD.Gs.size() - 2)];
  shared_ptr<Vertex> v = g->end;

  // fetch available time interval
  double
    t1 = g->start->position,
    t2 = v->G[1]->end->position;

  if (t2 - t1 <= 3*DBL_EPSILON) {
    cout << "DMC::shiftVertexPosition dt≈0 ->return" << endl;
    return;
  }

  // calculate exponent
  double
    c = (v->D[0]) ? -this->FD.phononEnergy(v->D[0]->q) : this->FD.phononEnergy(v->D[1]->q),
    dE = 0.5*v->G[0]->momentum.squaredNorm() - 0.5*v->G[1]->momentum.squaredNorm() - c,
    dtdE = (t2 - t1)*dE;

  // sample new t
  double t, r = this->Udouble(0, 1);
  if (-dtdE > 100) {
    // to avoid overflow due to exponential
    t = t2 - log(r)/dE;
  } else {
    t = t1 - log(1 - r*(1 - exp(-dtdE)))/dE;
  }

  double tOld = 0, oldVal = 0;
  if (this->debug) {
    tOld = v->position,
    oldVal = this->evaluateDiagram();
  }

  double boldContribution = 1;
  bool accepted = false;
  if (this->bold && this->boldIteration > 0) {
    // contribution from boldification
    boldContribution = exp(
                         this->additionalPhase(v->G[0]->p, t - v->G[0]->start->position)
                         + this->additionalPhase(v->G[1]->p, v->G[1]->end->position - t)
                         - this->additionalPhase(v->G[0])
                         - this->additionalPhase(v->G[1])
                       );

    // the rest of the acceptance ratio is unity
    if (boldContribution > this->Udouble(0, 1)) {
      this->FD.setVertexPosition(v, t);
      accepted = true;
    }

  } else {
    // is always accepted
    this->FD.setVertexPosition(v, t);
    accepted = true;
  }

  if ( ! isfinite(t)) {
    cout
      << "-------------------------" << endl
      << "DMC::shiftVertexPosition: nan encountered" << endl
      << "t=" << t << endl
      << "dtdE=" << dtdE << endl
      << "dE=" << dE << endl
      << "t1=" << t1 << endl
      << "t2=" << t2 << endl
      << "t2>t1=" << (t2>t1) << endl
      << "-------------------------" << endl;
  }

  if (this->debug) {
    double val = this->evaluateDiagram();

    double a;
    if (val == 0) {
      a = 0;
    } else if (oldVal == 0) {
      a = 1;
    } else {
      a = boldContribution * exp(dE*(tOld - t)) / (val/oldVal);
    }

    if (accepted) {
      this->checkAcceptanceRatio(a, "shiftVertexPosition");
    }

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::shiftVertexPosition " << endl
           << "a=" << a << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "dE=" << dE << endl
           << "exp=" << exp(dE*(t - tOld)) << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "shiftVertexPosition: " << a << endl;
    }
  }
}