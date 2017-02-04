#include "Phonon.h"

Phonon::Phonon (Vector3d _p, double theta, double phi) : Propagator(1, _p) {
  this->theta = theta;
  this->phi = phi;
}

void Phonon::print () {
  cout << "Phonon: " << this << endl
       << "\tmomentum: " << this->momentum.transpose() << endl
       << "\ttheta: " << this->theta << endl
       << "\tphi: " << this->phi << endl
       << "\tstart: " << this->start << endl
       << "\tend: " << this->end << endl;
}

void Phonon::setStart (Vertex *v) {
  // unlink propagator from previous vertex
  if (this->start) {
    start->setOutgoingD(NULL);
  }

  // unlink propagator linked to new vertex
  if (v->D[1] && v->D[1]->start) {
    v->D[1]->start = NULL;
  }

  // relink
  v->setOutgoingD(this);
  this->start = v;
}

void Phonon::setEnd (Vertex *v) {
  // unlink this propagator from previous vertex
  if (this->end) {
    this->end->setIngoingD(NULL);
  }

  // unlink propagator linked to new vertex
  if (v->D[0] && v->D[0]->end) {
    v->D[0]->end = NULL;
  }

  // relink
  v->setIngoingD(this);
  this->end = v;
}