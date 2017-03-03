#include "../DiagrammaticMonteCarloV2.h"

void DiagrammaticMonteCarloV2::changeDiagramLength (double param) {

  double
    l = 0.5*this->FD.end->G[0]->momentum.squaredNorm() - this->mu,
    r = this->Udouble(0, 1),
    tmin = this->FD.end->G[0]->start->position,
    dt = -log(1 - r + r*exp(-l*(this->maxLength - tmin)))/l;

  double oldVal, oldwInvt, wInvt;
  if (this->debug) {
    oldVal = this->FD();
    oldwInvt = exp(l*(this->FD.end->position - tmin))*(1 - exp(-l*this->maxLength - tmin))/l;
    wInvt = exp(l*dt)*(1 - exp(-l*this->maxLength - tmin))/l;
  }

  // is always accepted
  this->FD.setLength(tmin + dt);
  this->FD.end->setPosition(tmin + dt);

  if (this->debug) {
    double val = this->FD();

    double a;
    if (val == 0) {
      a = 0;
    } else if (oldwInvt == 0 || oldVal == 0) {
      a = 1;
    } else {
      a = val/oldVal * wInvt/oldwInvt;
    }

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::changeDiagramLength " << endl
           << "a=" << a << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "wInvt=" << wInvt << endl
           << "oldwInvt=" << oldwInvt << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "changeDiagramLength " << a << endl;
    }
  }
}