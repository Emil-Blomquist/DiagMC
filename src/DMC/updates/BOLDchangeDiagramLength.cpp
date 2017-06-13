#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::BOLDchangeDiagramLength (double param) {
  // select last electron propagator
  shared_ptr<Electron> g = this->FD.Gs.back();

  double
    l = this->lambdaOf(g),
    r = this->Udouble(0, 1),
    tmin = this->FD.end->G[0]->start->position;

  // phonon contribution when there are no external legs
  if (this->FD.Ds.size() > 0) {
    l += this->FD.phononEnergy(this->FD.end->D[0]->q);
  }

  double
    tofG = -log(1 - r + r*exp(-l*(this->maxLength - tmin)))/l,
    tofGold = g->end->position - g->start->position;

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
               + l*(tofG - tofGold)
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
      ratio = a * exp(-l*(tofG - tofGold)) / (val/oldVal);

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