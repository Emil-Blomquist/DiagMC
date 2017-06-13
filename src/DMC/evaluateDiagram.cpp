#include "DiagrammaticMonteCarlo.h"

double DiagrammaticMonteCarlo::evaluateDiagram () {
  double val = 1;

  if (this->bold && this->boldIteration > 0) {
    for (auto g : this->FD.Gs) {
      val *= (*g)(this->mu, this->dEOf(g));
    }
  } else {
    for (auto g : this->FD.Gs) {
      val *= (*g)(this->mu);
    }
  }

  for (auto d : this->FD.Ds) {
    val *= (*d)(this->alpha, this->FD.phononEnergy(d->q));
  }

  return val;
}