#include "Electron.h"

Electron::Electron (Vector3d P) : Propagator(0, P) {}

void Electron::print () {
  cout << "Electron: " << this->shared_from_this() << endl
       << "\tmomentum: " << this->momentum.transpose() << endl
       << "\tstart: " << this->start << endl
       << "\tend: " << this->end << endl;
}

void Electron::setStart (shared_ptr<Vertex> v) {
  // unlink propagator from previous vertex
  if (this->start) {
    start->setG(1, NULL);
  }

  // unlink propagator linked to new vertex
  if (v && v->G[1] && v->G[1]->start) {
    v->G[1]->start.reset();
  }

  // relink
  if (v) v->setG(1, this->shared_from_this());
  this->start = v;
}

void Electron::setEnd (shared_ptr<Vertex> v) {
  // unlink this propagator from previous vertex
  if (this->end) {
    this->end->setG(0, NULL);
  }

  // unlink propagator linked to new vertex
  if (v && v->G[0] && v->G[0]->end) {
    v->G[0]->end.reset();
  }

  // relink
  if (v) v->setG(0, this->shared_from_this());
  this->end = v;
}