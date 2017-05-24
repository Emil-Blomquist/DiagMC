#include "Electron.h"

Electron::Electron (Vector3d P) {
  this->momentum = P;
  this->p = P.norm();
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

void Electron::addMomentum (Vector3d P) {
  this->momentum += P;

  this->p = this->momentum.norm();

  if ( ! isfinite(P[0]) || ! isfinite(P[1]) || ! isfinite(P[2])) {
    cout << "Phonon::addMomentum P=" << P.transpose() << endl;
  }
}

double Electron::operator() (double mu, double dE) {
  double
    E = 0.5*pow(this->p, 2.0) + dE - mu,
    t = this->end->position - this->start->position;

  return exp(-E*t);
}