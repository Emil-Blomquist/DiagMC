#include "../DiagrammaticMonteCarloV1.h"

double DiagrammaticMonteCarloV1::lowerOrder (double param) {
  // requirement to lower: must be at least of order 1
  if (this->FD.Ds.size() == 0) {
    return 0;
  }

  // old configuration value
  double oldVal = this->FD();

  // pick a phonon line
  int i = this->Uint(0, this->FD.Ds.size() - 1);
  shared_ptr<Phonon> d = this->FD.Ds[i];
  double Winvd = this->FD.Ds.size();

  // vertices
  shared_ptr<Vertex>
    v1 = d->start,
    v2 = d->end;

  // calculate vertex probability
  int i1 = 0, i2 = 0;
  for (int i = 0; i != this->FD.Gs.size(); ++i) {
    if (this->FD.Gs[i]->end == v1) {
      i1 = i;
      break;
    }
  }
  for (int i = i1 + 1; i != this->FD.Gs.size(); ++i) {
    if (this->FD.Gs[i]->end == v2) {
      i2 = i;
      break;
    }
  }
  double WinvVertex = (this->FD.Gs.size() - 2) * (this->FD.Gs.size() - 2 - i1);

  // inverse probabilities for internal parameters
  // must be the same probabilities as for raising order!
  double
    wInvt1 = (v1->G[1]->end == v2 ? v1->G[1]->end->G[1]->end->position : v1->G[1]->end->position) - v1->G[0]->start->position,

    l = 4,
    Dt2 = v2->G[1]->end->position - v2->G[0]->start->position,
    dt2 = v2->position - v2->G[0]->start->position,
    wInvt2 = exp(l*dt2)*(1 - exp(-l*Dt2))/l,

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
    a = val/oldVal * Winvd/(WinvVertex*wInvt1*wInvt2*wInvQ);
  }

  if (this->debug) {
    if (val == 0 || wInvQ == 0) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::lowerOrder " << this->FD.Ds.size() << " " << a << endl
           << "--------------------------------------------------------------------" << endl;
    } else {
      if (this->loud) { cout << "lowerOrder: "  << this->FD.Ds.size() << " " << a << endl; }
    }
  }
  
  // these should be unlinked if update is accepted
  this->phonon2beRemoved = d;
  this->vertices2beRemoved[0] = v1;
  this->vertices2beRemoved[1] = v2;
  this->electrons2beRemoved[0] = v1->G[1];
  this->electrons2beRemoved[1] = v2->G[1];
  
  return a;
}