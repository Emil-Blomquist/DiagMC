#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::BOLDshiftVertexPosition (double param) {
  cout << "BOLDshiftVertexPosition: Do use the ordinary one instead. This one has pretty bad statistics" << endl;


  // requirement to lower: must be at least of order 1
  if ( this->FD.Ds.size() == 0 || this->FD.Ds.size() == 1) {
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
    cout << "DMC::BOLDshiftVertexPosition dt≈0 ->return" << endl;
    return;
  }

  // calculate exponent
  double
    c = (v->D[0]) ? -this->FD.phononEnergy(v->D[0]->q) : this->FD.phononEnergy(v->D[1]->q),
    dE = 0.5*v->G[0]->momentum.squaredNorm() - 0.5*v->G[1]->momentum.squaredNorm() - c,
    dtdE = (t2 - t1)*dE;

  // sample new t
  double t = this->Udouble(t1, t2);

  double tOld = 0, oldVal = 0;
  if (this->debug) {
    tOld = v->position,
    oldVal = this->evaluateDiagram();
  }


  // contribution from boldification
  double boldContribution = this->additionalPhase(v->G[0]->p, t - v->G[0]->start->position)
                          + this->additionalPhase(v->G[1]->p, v->G[1]->end->position - t)
                          - this->additionalPhase(v->G[0])
                          - this->additionalPhase(v->G[1]);


  // acceptance ratio
  double a = exp(boldContribution + dE*(tOld - t));

  bool accepted = false;
  if (a > this->Udouble(0, 1)) {
    this->FD.setVertexPosition(v, t);
    accepted = true;
  }

  if ( ! isfinite(t)) {
    cout
      << "-------------------------" << endl
      << "DMC::BOLDshiftVertexPosition: nan encountered" << endl
      << "t=" << t << endl
      << "dtdE=" << dtdE << endl
      << "dE=" << dE << endl
      << "t1=" << t1 << endl
      << "t2=" << t2 << endl
      << "t2>t1=" << (t2>t1) << endl
      << "-------------------------" << endl;
  }

  if (this->debug) {
    double
      val = this->evaluateDiagram(),
      ratio = a / (val/oldVal);


    if (accepted) {
      this->checkAcceptanceRatio(ratio, "BOLDshiftVertexPosition");
    }

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::BOLDshiftVertexPosition " << endl
           << "ratio=" << ratio << endl
           << "a=" << a << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "dE=" << dE << endl
           << "exp=" << exp(dE*(t - tOld)) << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "BOLDshiftVertexPosition: " << a << endl;
    }
  }
}