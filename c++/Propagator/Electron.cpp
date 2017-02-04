#include "Electron.h"

Electron::Electron (Vector3d _p) : Propagator(0, _p) {}

void Electron::print () {
  cout << "Electron: " << this << endl
       << "\tmomentum: " << this->momentum.transpose() << endl
       << "\tstart: " << this->start << endl
       << "\tend: " << this->end << endl;
}

void Electron::setStart (Vertex *v) {
  // unlink propagator from previous vertex
  if (this->start) {
    start->setOutgoingG(NULL);
  }

  // unlink propagator linked to new vertex
  if (v->G[1] && v->G[1]->start) {
    v->G[1]->start = NULL;
  }

  // relink
  v->setOutgoingG(this);
  this->start = v;
}

void Electron::setEnd (Vertex *v) {
  // unlink this propagator from previous vertex
  if (this->end) {
    this->end->setIngoingG(NULL);
  }

  // unlink propagator linked to new vertex
  if (v->G[0] && v->G[0]->end) {
    v->G[0]->end = NULL;
  }

  // relink
  v->setIngoingG(this);
  this->end = v;
}