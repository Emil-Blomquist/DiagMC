#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::BOLDchangeDiagramLength (double param) {
  // select last electron propagator
  shared_ptr<Electron> g = this->FD.Gs.back();

  double
    tmin = g->start->position,
    tofG = this->Udouble(0, this->maxLength - tmin);

  double oldVal = 0;
  if (this->debug) {
    oldVal = this->evaluateDiagram();
  }

  // phonon energy contribution
  double phononEnergy = 0;
  if (this->FD.Ds.size() > 0) {
    phononEnergy = this->FD.phononEnergy(this->FD.end->D[0]->q);
  }

  double a = exp(
               (this->mu - 0.5*pow(g->p, 2.0) - phononEnergy)*(tofG - (this->FD.length - tmin))
               + this->additionalPhase(g->p, tofG)
               - this->additionalPhase(g)
             );

  bool accepted = false;
  if (a > this->Udouble(0, 1)) {
    this->FD.setLength(tmin + tofG);
    this->FD.end->setPosition(tmin + tofG);
    accepted = true;
  }

  if (this->debug) {
    double
      val = this->evaluateDiagram(),
      ratio = a / (val/oldVal);

    if (accepted) {
      this->checkAcceptanceRatio(ratio, "BOLDchangeDiagramLength");
    }

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::BOLDchangeDiagramLength " << endl
           << "a=" << a << endl 
           << "ratio=" << ratio << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "BOLDchangeDiagramLength: " << a << endl;
    }
  }
}