#include "Phonon.h"

Phonon::Phonon (Vector3d Q, double q, double theta, double phi) {
  this->momentum = Q;
  this->q = q;
  this->theta = theta;
  this->phi = phi;
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

void Phonon::setMomentum (Vector3d Q, double q) {
  this->momentum = Q;
  this->q = q;

  if ( ! isfinite(q)) {
    cout << "Phonon::setMomentum Q=" << Q.transpose() << endl;
  }
}

void Phonon::setTheta (double theta) {
  this->theta = theta;

  if ( ! isfinite(theta)) {
    cout << "Phonon::setTheta theta=" << theta << endl;
  }
}

void Phonon::setPhi (double phi) {
  this->phi = phi;

  if ( ! isfinite(phi)) {
    cout << "Phonon::setPhi phi=" << phi << endl;
  }
}

double Phonon::operator() (double alpha, double omega) {
  double t = this->end->position - this->start->position;

  return this->value(omega, t, this->theta, alpha);
}

double Phonon::value (double omega, double t, double theta, double alpha) {
  return alpha/(sqrt(8)*M_PI*M_PI) * sin(theta) * exp(-omega*t);
}