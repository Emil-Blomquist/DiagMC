#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarloV2::shiftVertexPosition (double param) {
  // requirement to lower: must be at least of order 1
  if (this->FD.Ds.size() == 0) {
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
    c = (v->D[0]) ? -1 : 1,
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
    oldVal = this->FD();
  }

  // is always accepted
  this->FD.setVertexPosition(v, t);

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
    double val = this->FD();

    double a;
    if (val == 0) {
      a = 0;
    } else if (oldVal == 0) {
      a = 1;
    } else {
      a = val/oldVal * exp(dE*(t - tOld));
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
      cout << "shiftVertexPosition " << a << endl;
    }
  }
}