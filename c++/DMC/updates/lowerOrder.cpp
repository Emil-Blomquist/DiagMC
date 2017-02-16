#include "../DiagrammaticMonteCarlo.h"

double DiagrammaticMonteCarlo::lowerOrder () {
  // requirement to lower: must be at least of order 1
  if (this->FD.Ds.size() == 0) {
    return 0;
  }

  // first encountered phonon line
  shared_ptr<Phonon> d = this->FD.Gs[0]->end->D[1];

  // requirement to lower: first and third vertex must connect to the same phonon
  if (d != this->FD.Gs[2]->end->D[0]) {
    return 0;
  }

  // old configuration value
  double oldVal = this->FD();

  // vertices
  shared_ptr<Vertex>
    v1 = d->start,
    v2 = d->end;

  // inverse probabilities for internal parameters
  // must be the same probabilities as for raising order!
  double wInvt1, wInvt2;
  if (this->FD.Ds.size() == 1) {
    wInvt1 = this->FD.length;
    wInvt2 = d->end->G[1]->end->position - d->end->G[0]->start->position;
  } else {
    wInvt1 = d->start->G[1]->end->position - d->start->G[0]->start->position;
    wInvt2 = d->end->G[1]->end->position - d->end->G[0]->start->position;
  }

  double
    std = 1/sqrt(v2->position - v1->position),
    wInvQ = 2*pow(M_PI, 2.0) * sqrt(0.5*M_PI*pow(std, 2.0)) * exp(0.5*pow(d->momentum.norm()/std, 2.0));

  // remove phonon
  this->FD.removeInternalPhonon(d);

  // remove vertices
  this->FD.removeVertex(v1);
  this->FD.removeVertex(v2);

  // current configuration value
  double val = this->FD();

  // acceptance ration
  double a;
  if (val == 0 || wInvQ == 0) {
    a = 1;
  } else {
    a = val/oldVal / (wInvt1 * wInvt2 * wInvQ);
  }

  if (this->debug) {
    if (val == 0 || wInvQ == 0) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::lowerOrder " << a << endl
           << "--------------------------------------------------------------------" << endl;
    } else {
      cout << "lowerOrder: " << a << endl;
    }
  }
  
  return a;
}


// # get current diagram value
// diag = self.FD()

// # acceptance ratio
// # if diagOld == 0 or wInvQ == 0:
// #   R = 1
// # else:
// R = diag/diagOld / (wInvt1 * wInvt2 * wInvQ)

// return R