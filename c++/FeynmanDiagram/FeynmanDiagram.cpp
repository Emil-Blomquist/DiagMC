#include "FeynmanDiagram.h"

FeynmanDiagram::FeynmanDiagram (Vector3d _externalMomentum, double _length, double _couplingConstant, double _chemicalPotential) {
  this->externalMomentum = _externalMomentum;
  this->length = _length;
  this->couplingConstant = _couplingConstant;
  this->chemicalPotential = _chemicalPotential;

  // initiate bare propagator
  Vertex v1(0), v2(length);
  this->start = &v1; this->end = &v2;


  Electron g(externalMomentum);
  g.setStart(start);
  g.setEnd(end);
  this->Gs.push_back(&g);
}