#include "Propagator.h"

Propagator::Propagator (int type, Vector3d momentum) {
  this->type = type;
  this->momentum = momentum;
}

void Propagator::setMomentum (Vector3d P) {
  this->momentum = P;
}

void Propagator::addMomentum (Vector3d P) {
  this->momentum += P;
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