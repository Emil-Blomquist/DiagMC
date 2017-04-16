#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::changeDiagramLength (double param) {

  double
    l = 0.5*this->FD.end->G[0]->momentum.squaredNorm() - this->mu,
    r = this->Udouble(0, 1),
    tmin = this->FD.end->G[0]->start->position;

  // phonon contribution
  if ( ! this->externalLegs && this->FD.Ds.size() > 0) {
    l += this->FD.phononEnergy(this->FD.end->D[0]->q);
  }

  double dt = -log(1 - r + r*exp(-l*(this->maxLength - tmin)))/l;

  double oldVal = 0, oldwInvt = 0, wInvt = 0;
  if (this->debug) {
    oldVal = this->FD();
    oldwInvt = exp(l*(this->FD.end->position - tmin))*(1 - exp(-l*(this->maxLength - tmin)))/l;
    wInvt = exp(l*dt)*(1 - exp(-l*(this->maxLength - tmin)))/l;
  }

  // is always accepted
  this->FD.setLength(tmin + dt);
  this->FD.end->setPosition(tmin + dt);

  if ( ! isfinite(dt)) {
    cout
      << "-------------------------" << endl
      << "DMC::changeDiagramLength: nan encountered" << endl
      << "dt=" << dt << endl
      << "tmin=" << tmin << endl
      << "r=" << r << endl
      << "l=" << l << endl
      << "-------------------------" << endl;
  }

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

    this->checkAcceptanceRatio(a, "changeDiagramLength");

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
      cout << "changeDiagramLength: " << a << endl;
    }
  }
}