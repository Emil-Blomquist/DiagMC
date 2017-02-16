#include "Phonon.h"

Phonon::Phonon (Vector3d _p, double theta, double phi) : Propagator(1, _p) {
  this->theta = theta;
  this->phi = phi;
}

void Phonon::print () {
  cout << "Phonon: " << this->shared_from_this() << endl
       << "\tmomentum: " << this->momentum.transpose() << endl
       << "\ttheta: " << this->theta << endl
       << "\tphi: " << this->phi << endl
       << "\tstart: " << this->start << endl
       << "\tend: " << this->end << endl;
}

void Phonon::setStart (shared_ptr<Vertex> v) {
  // unlink propagator from previous vertex
  if (this->start) {
    this->start->setD(1, NULL);
  }

  // unlink propagator linked to new vertex
  if (v && v->D[1] && v->D[1]->start) {
    v->D[1]->start.reset();
  }

  // relink
  if (v) v->setD(1, this->shared_from_this());
  this->start = v;
}

void Phonon::setEnd (shared_ptr<Vertex> v) {
  // unlink this propagator from previous vertex
  if (this->end) {
    this->end->setD(0, NULL);
  }

  // unlink propagator linked to new vertex
  if (v && v->D[0] && v->D[0]->end) {
    v->D[0]->end.reset();
  }

  // relink
  if (v) v->setD(0, this->shared_from_this());
  this->end = v;
}

void Phonon::setTheta (double theta) {
  this->theta = theta;
}

void Phonon::setPhi (double phi) {
  this->phi = phi;
}

double Phonon::operator() (double alpha) {
  double E, t;
  E = 1;
  t = this->end->position - this->start->position;

  return alpha/(sqrt(8)*M_PI*M_PI) * sin(this->theta) * exp(-E*t);
}

void Phonon::save () {
  Propagator::save();
  this->savedTheta = this->theta;
  this->savedPhi = this->phi;
}

void Phonon::revert () {
  Propagator::revert();
  this->theta = this->savedTheta;
  this->phi = this->savedPhi;
}