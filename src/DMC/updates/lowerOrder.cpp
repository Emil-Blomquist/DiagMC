#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::lowerOrder (double param) {
  // requirement to lower: must be at least of order 1
  if (this->FD.Ds.size() == 0) {
    return;
  }

  // old configuration value
  double oldVal = 0;
  if (this->debug) {
    oldVal = this->FD();
  }

  // pick a phonon line
  int phononIndex = this->Uint(0, this->FD.Ds.size() - 1);
  shared_ptr<Phonon> d = this->FD.Ds[phononIndex];
  double Winvd = this->FD.Ds.size();

  // vertices
  shared_ptr<Vertex>
    v1 = d->start,
    v2 = d->end;

  // calculate vertex probability
  double WinvVertex = this->FD.Gs.size() - 2;

  // inverse probabilities for internal parameters
  // must be the same probabilities as for raising order!
  double
    wInvt1 = (v1->G[1]->end == v2 ? v2->G[1]->end->position : v1->G[1]->end->position) - v1->G[0]->start->position,

    t2low = v1->position,
    t2up = this->FD.length,
    l = 0.01,
    dt2 = v2->position - t2low,
    wInvt2 = exp(l*dt2)*(1 - exp(-l*(t2up - t2low)))/l,

    std = 1/sqrt(v2->position - v1->position),
    wInvQ = 2*pow(M_PI, 2.0) * sqrt(0.5*M_PI*pow(std, 2.0)) * exp(0.5*pow(d->momentum.norm()/std, 2.0));

  // to calculate the acceptance ratio
  Vector3d
    Q = d->momentum,
    P0 = this->calculateP0(d);

  double
    sinTheta = sin(d->theta),
    alpha = this->alpha,
    q2 = Q.squaredNorm(),
    dt = v2->position - v1->position,
    exponential = exp(dt*(1 + 0.5*q2 - Q.dot(P0)));

  double a;
  if (wInvQ == 0 || sinTheta == 0 || wInvt2 == 0 || ! isfinite(exponential)) {
    a = 1;
  } else {
    a = exponential * sqrt(8)*M_PI*M_PI/(alpha*sinTheta) * Winvd/(WinvVertex*wInvt1*wInvt2*wInvQ);
  }

  // accept or reject update
  bool accepted = false;
  if (a > this->Udouble(0, 1)) {
    // remove phonon
    this->FD.removeInternalPhonon(phononIndex);

    // remove vertices
    this->FD.removeVertex(v1);
    this->FD.removeVertex(v2);

    accepted = true;
  }

  if (this->debug) {
    double val = this->FD();

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::lowerOrder " << endl
           << "accepted=" << accepted << endl
           << "a=" << a << endl
           << "a_diag" << val/oldVal * Winvd/(WinvVertex*wInvt1*wInvt2*wInvQ) << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "wInvQ=" << wInvQ << endl
           << "wInvt2=" << wInvt2 << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "lowerOrder: "  << this->FD.Ds.size() << " " << a << endl;
    }
  }
}










// void DiagrammaticMonteCarlo::lowerOrder (double param) {
//   // requirement to lower: must be at least of order 1
//   if (this->FD.Ds.size() == 0) {
//     return;
//   }

//   // old configuration value
//   double oldVal;
//   if (this->debug) {
//     oldVal = this->FD();
//   }

//   // pick a phonon line
//   int phononIndex = this->Uint(0, this->FD.Ds.size() - 1);
//   shared_ptr<Phonon> d = this->FD.Ds[phononIndex];
//   double Winvd = this->FD.Ds.size();

//   // vertices
//   shared_ptr<Vertex>
//     v1 = d->start,
//     v2 = d->end;

//   // calculate vertex probability
//   int i1 = 0;
//   for (unsigned int i = 0; i != this->FD.Gs.size(); ++i) {
//     if (this->FD.Gs[i]->end == v1) {
//       i1 = i;
//       break;
//     }
//   }
//   double WinvVertex = (this->FD.Gs.size() - 2) * (this->FD.Gs.size() - 2 - i1);

//   // inverse probabilities for internal parameters
//   // must be the same probabilities as for raising order!
//   double
//     wInvt1 = (v1->G[1]->end == v2 ? v1->G[1]->end->G[1]->end->position : v1->G[1]->end->position) - v1->G[0]->start->position,

//     l = 4,
//     Dt2 = v2->G[1]->end->position - v2->G[0]->start->position,
//     dt2 = v2->position - v2->G[0]->start->position,
//     wInvt2 = exp(l*dt2)*(1 - exp(-l*Dt2))/l,

//     std = 1/sqrt(v2->position - v1->position),
//     wInvQ = 2*pow(M_PI, 2.0) * sqrt(0.5*M_PI*pow(std, 2.0)) * exp(0.5*pow(d->momentum.norm()/std, 2.0));

//   // to calculate the acceptance ratio
//   Vector3d
//     Q = d->momentum,
//     P0 = this->calculateP0(d);

//   double
//     sinTheta = sin(d->theta),
//     alpha = this->alpha,
//     q2 = Q.squaredNorm(),
//     dt = v2->position - v1->position,
//     exponential = exp(dt*(1 + 0.5*q2 - Q.dot(P0)));

//   double a;
//   if (wInvQ == 0 || sinTheta == 0 || wInvt2 == 0 || ! isfinite(exponential)) {
//     a = 1;
//   } else {
//     a = exponential * sqrt(8)*M_PI*M_PI/(alpha*sinTheta) * Winvd/(WinvVertex*wInvt1*wInvt2*wInvQ);
//   }

//   // accept or reject update
//   bool accepted = false;
//   if (a > this->Udouble(0, 1)) {
//     // remove phonon
//     this->FD.removeInternalPhonon(phononIndex);

//     // remove vertices
//     this->FD.removeVertex(v1);
//     this->FD.removeVertex(v2);

//     accepted = true;
//   }

//   if (this->debug) {
//     double val = this->FD();

//     if (a < 0 || ! isfinite(a)) {
//       cout << "--------------------------------------------------------------------" << endl
//            << "overflow at DMC::lowerOrder " << endl
//            << "accepted=" << accepted << endl
//            << "a=" << a << endl
//            << "a_diag" << val/oldVal * Winvd/(WinvVertex*wInvt1*wInvt2*wInvQ) << endl
//            << "order=" << this->FD.Ds.size() << endl
//            << "val=" << val << endl
//            << "oldVal=" << oldVal << endl
//            << "wInvQ=" << wInvQ << endl
//            << "wInvt2=" << wInvt2 << endl
//            << "--------------------------------------------------------------------" << endl;
//     } else if (this->loud) {
//       cout << "lowerOrder: "  << this->FD.Ds.size() << " " << a << endl;
//     }
//   }
// }