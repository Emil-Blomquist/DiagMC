#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::changeDiagramLengthComplex (double param) {
  // pick a electron propagator on random
  shared_ptr<Electron> g = this->FD.Gs[this->Uint(0, this->FD.Gs.size() - 1)];

  double phononContribution = 0;

  if (
    (g->start->G[0] && g->start->G[0]->start->position == g->start->position) ||
    (g->end->G[1] && g->end->G[1]->end->position == g->end->position) ||
    g->start->position == g->end->position
  ) {
    // a neighbour has equal times

    cout << "DMC::changeDiagramLengthComplex: equal time fix" << endl;

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
          phononContribution += 1;
        }
      }
    }
  } else {
    // neighbours have separate times
    for (shared_ptr<Phonon> d : this->FD.Ds) {
      if (d->start->position <= g->start->position && d->end->position >= g->end->position) {
        phononContribution += 1;
      }
    }
  }

  double
    l = 0.5*g->momentum.squaredNorm() + phononContribution - this->mu,
    r = this->Udouble(0, 1),
    lowerBound = g->start->position,
    upperBound = this->maxLength - (this->FD.length - g->end->position),
    gLength = -log(1 - r + r*exp(-l*(upperBound - lowerBound)))/l,
    dt = gLength - (g->end->position - g->start->position);

  double oldVal = 0, oldwInvt = 0, wInvt = 0;
  if (this->debug) {
    oldVal = this->FD();
    oldwInvt = exp(l*(g->end->position - g->start->position))*(1 - exp(-l*(upperBound - lowerBound)))/l;
    wInvt = exp(l*gLength)*(1 - exp(-l*(upperBound - lowerBound)))/l;
  }

  // is always accepted
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

  if ( ! isfinite(dt)) {
    cout
      << "-------------------------" << endl
      << "DMC::changeDiagramLengthComplex: nan encountered" << endl
      << "gLength=" << gLength << endl
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

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::changeDiagramLengthComplex " << endl
           << "a=" << a << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "wInvt=" << wInvt << endl
           << "oldwInvt=" << oldwInvt << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "changeDiagramLengthComplex " << a << endl;
    }
  } 
}