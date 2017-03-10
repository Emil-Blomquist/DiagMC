// #include "../DiagrammaticMonteCarloV1.h"

// double DiagrammaticMonteCarloV1::raiseOrder (double param) {
//   // requirement to raise: cannot be of max order
//   if (this->FD.Ds.size() == this->maxOrder) {
//     return 0;
//   }

//   // old configuration value
//   double oldVal = this->FD();

//   // select electron lines to split
//   int
//     i1 = this->Uint(0, this->FD.Gs.size() - 1),
//     i2 = this->Uint(i1 + 1, this->FD.Gs.size());
//   double wInvVertex = this->FD.Gs.size() * (this->FD.Gs.size() - i1);

//   // split first electron line
//   shared_ptr<Electron> g1 = this->FD.Gs[i1];
//   double
//     t1 = this->Udouble(g1->start->position, g1->end->position),
//     wInvt1 = g1->end->position - g1->start->position;
//   shared_ptr<Vertex> v1 = this->FD.insertVertex(i1, t1 - g1->start->position);

//   // split second electron line
//   shared_ptr<Electron> g2 = this->FD.Gs[i2];
//   double
//     l = 4,
//     Dt2 = g2->end->position - g2->start->position,
//     r = this->Udouble(0, 1),
//     dt2 = -log(1 - r + r*exp(-l*Dt2))/l,
//     wInvt2 = exp(l*dt2)*(1 - exp(-l*Dt2))/l;
//   shared_ptr<Vertex> v2 = this->FD.insertVertex(i2, dt2);

//   // generate phonon momentum
//   double std = 1/sqrt(g2->start->position + dt2 - t1);
//   normal_distribution<double> normal(0.0, std);

//   double
//     q = abs(normal(this->mt)),
//     theta = this->Udouble(0, M_PI),
//     phi = this->Udouble(0, 2*M_PI),
//     wInvQ = 2*pow(M_PI, 2.0) * sqrt(0.5*M_PI*pow(std, 2.0)) * exp(0.5*pow(q/std, 2.0));

//   // phonon
//   double wInvd = this->FD.Ds.size() + 1;

//   Vector3d Q(q*sin(theta)*cos(phi), q*sin(theta)*sin(phi), q*cos(theta));

//   // add phonon
//   shared_ptr<Phonon> d = this->FD.addInternalPhonon(v1, v2, Q, theta, phi);

//   // current confiduration value
//   double val = this->FD();

//   // acceptance ration
//   double a = val/oldVal * (wInvVertex*wInvt1*wInvt2*wInvQ)/wInvd;
//   // double a = (wInvVertex*wInvt1)/wInvd;

//   if (this->debug) {
//     if (this->loud) { cout << "raiseOrder: " << a << endl; }
//   }

//   // these should be unlinked if update is rejected
//   this->phonon2beRemoved = d;
//   this->vertices2beRemoved[0] = v1;
//   this->vertices2beRemoved[1] = v2;
//   this->electrons2beRemoved[0] = v1->G[1];
//   this->electrons2beRemoved[1] = v2->G[1];

//   return a;
// }