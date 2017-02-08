#include "FeynmanDiagram.h"

FeynmanDiagram::FeynmanDiagram (
  Vector3d _externalMomentum,
  double _length,
  double _couplingConstant,
  double _chemicalPotential
) {
  this->externalMomentum = _externalMomentum;
  this->length = _length;
  this->couplingConstant = _couplingConstant;
  this->chemicalPotential = _chemicalPotential;

  // create vertices
  this->start.reset(new Vertex(0));
  this->end.reset(new Vertex(length));
  
  // create electronic propagator
  this->Gs.emplace_back(new Electron(externalMomentum));

  // link propagator to vertices
  this->Gs[0]->setStart(this->start);
  this->Gs[0]->setEnd(this->end);
}

void FeynmanDiagram::print () {
  this->start->print();
  this->end->print();

  for (auto g : this->Gs) {
    g->print();
  }

  for (auto d : this->Ds) {
    d->print();
  }
}