#include "Propagator.h"

Propagator::Propagator (int type, Vector3d momentum) {
  this->type = type;
  this->momentum = momentum;
}

void Propagator::setMomentum (Vector3d P) {
  this->momentum = P;

  if ( ! isfinite(P[0]) || ! isfinite(P[1]) || ! isfinite(P[2])) {
    cout << "Phonon::setMomentum P=" << P.transpose() << endl;
  }
}

void Propagator::addMomentum (Vector3d P) {
  this->momentum += P;

  if ( ! isfinite(P[0]) || ! isfinite(P[1]) || ! isfinite(P[2])) {
    cout << "Phonon::addMomentum P=" << P.transpose() << endl;
  }
}

void Propagator::save () {
  this->savedMomentum = this->momentum;
  this->savedStart = this->start;
  this->savedEnd = this->end;
}

void Propagator::revert () {
  this->momentum = this->savedMomentum;
  this->start = this->savedStart;
  this->end = this->savedEnd;
}

void Propagator::unlink () {
  this->start = NULL;
  this->end = NULL;
  this->savedStart = NULL;
  this->savedEnd = NULL;
}