#include "../DiagrammaticMonteCarlo.h"

double DiagrammaticMonteCarlo::raiseOrder () {
  // requirement to raise: cannot be of max order
  if (this->FD.Ds.size() == this->maxOrder) {
    return 0;
  }

  // old configuration value
  double oldVal = this->FD();

  // select electron lines to split
  int i1, i2;
  if (this->FD.Ds.size() == 0) {
    i1 = 0; i2 = 1;
  } else {
    i1 = 0; i2 = 2;
  }

  // split first electron line
  shared_ptr<Electron> g1 = this->FD.Gs[i1];
  double
    t1 = this->Udouble(g1->start->position, g1->end->position),
    wInvt1 = g1->end->position - g1->start->position;
  shared_ptr<Vertex> v1 = this->FD.insertVertex(i1, t1 - g1->start->position);

  // split second electron line
  shared_ptr<Electron> g2 = this->FD.Gs[i2];
  double
    t2 = this->Udouble(g2->start->position, g2->end->position),
    wInvt2 = g2->end->position - g2->start->position;
  shared_ptr<Vertex> v2 = this->FD.insertVertex(i2, t2 - g2->start->position);

  // generate phonon momentum
  double std = 1/sqrt(t2 - t1);
  normal_distribution<double> normal(0.0, std);

  double
    q = abs(normal(this->mt)),
    theta = this->Udouble(0, M_PI),
    phi = this->Udouble(0, 2*M_PI),
    wInvQ = 2*pow(M_PI, 2.0) * sqrt(0.5*M_PI*pow(std, 2.0)) * exp(0.5*pow(q/std, 2.0));

  Vector3d Q(q*sin(theta)*cos(phi), q*sin(theta)*sin(phi), q*cos(theta));

  // add phonon
  shared_ptr<Phonon> d = this->FD.addInternalPhonon(v1, v2, Q, theta, phi);

  // current confiduration value
  double val = this->FD();

  // acceptance ration
  double a = val/oldVal * (wInvt1 * wInvt2 * wInvQ);

  if (this->debug) {
    if (this->loud) { cout << "raiseOrder: " << a << endl; }
  }

  // these should be unlinked if update is rejected
  this->phonon2beRemoved = d;
  this->vertices2beRemoved[0] = v1;
  this->vertices2beRemoved[1] = v2;
  this->electrons2beRemoved[0] = v1->G[1];
  this->electrons2beRemoved[1] = v2->G[1];

  return a;
}
