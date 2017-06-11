#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::BOLDchangeDiagramLengthComplex (double param) {
  // pick a electron propagator on random
  shared_ptr<Electron> g = this->FD.Gs[this->Uint(0, this->FD.Gs.size() - 1)];

  double phononEnergies = 0;

  if (
    (g->start->G[0] && g->start->G[0]->start->position == g->start->position) ||
    (g->end->G[1] && g->end->G[1]->end->position == g->end->position) ||
    g->start->position == g->end->position
  ) {
    // a neighbour has equal times

    cout << "DMC::BOLDchangeDiagramLengthComplex: equal time fix" << endl;

    for (shared_ptr<Phonon> d : this->FD.Ds) {
      shared_ptr<Vertex> v = this->FD.start;

      bool foundFirstVertex = false;
      while (v != g->end) {
        if (v == d->start) {
          foundFirstVertex = true;
          break;
        }

        v = v->G[1]->end;
      }

      if (foundFirstVertex) {

        v = g->end;

        bool foundLastVertex = false;
        while (v != this->FD.end) {
          if (v == d->end) {
            foundLastVertex = true;
            break;
          }

          v = v->G[1]->end;
        }

        if (foundLastVertex) {
          phononEnergies += this->FD.phononEnergy(d->q);
        }
      }
    }
  } else {
    // neighbours have separate times
    for (shared_ptr<Phonon> d : this->FD.Ds) {
      if (d->start->position <= g->start->position && d->end->position >= g->end->position) {
        phononEnergies += this->FD.phononEnergy(d->q);
      }
    }
  }

   double
    l = this->lambdaOf(g) + phononEnergies,
    r = this->Udouble(0, 1),
    lowerBound = g->start->position,
    upperBound = this->maxLength - (this->FD.length - g->end->position),
    tofG = -log(1 - r + r*exp(-l*(upperBound - lowerBound)))/l,
    tofGold = g->end->position - g->start->position,
    dt = tofG - (g->end->position - g->start->position);

  double oldVal = 0;
  if (this->debug) {
    oldVal = this->evaluateDiagram();
  }

  double a = exp(
               (this->mu - 0.5*pow(g->p, 2.0) - phononEnergies)*(tofG - (g->end->position - g->start->position))
               + this->additionalPhase(g->p, tofG)
               - this->additionalPhase(g)
               + l*(tofG - tofGold)
             );

  bool accepted = false;
  if (a > this->Udouble(0, 1)) {
    this->FD.setLength(this->FD.length + dt);

    if (dt > 0) {
      // increasing time -> start from back
      shared_ptr<Vertex> v = this->FD.end;

      while (v != g->start) {
        v->setPosition(v->position + dt);

        v = v->G[0]->start;
      }
    } else {
      // decreasing time -> start from front
      shared_ptr<Vertex> v = g->start;

      do {
        v = v->G[1]->end;

        v->setPosition(v->position + dt);
      } while (v->G[1]);
    }

    accepted = true;
  }

  if (this->debug) {
    double
      val = this->evaluateDiagram(),
      ratio = a * exp(-l*(tofG - tofGold)) / (val/oldVal);

    if (accepted) {
      this->checkAcceptanceRatio(ratio, "BOLDchangeDiagramLengthComplex");
    }

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::BOLDchangeDiagramLengthComplex " << endl
           << "a=" << a << endl 
           << "ratio=" << ratio << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "BOLDchangeDiagramLengthComplex: " << a << endl;
    }
  }
}