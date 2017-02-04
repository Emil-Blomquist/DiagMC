#include "Propagator.h"

Propagator::Propagator (int type, Vector3d momentum) {
  this->type = type;
  this->momentum = momentum;

  this->start = NULL;
  this->end = NULL;
}

void Propagator::setMomentum (Vector3d P) {
 this->momentum = P;
}

void Propagator::addMomentum (Vector3d P) {
 this->momentum += P;
}