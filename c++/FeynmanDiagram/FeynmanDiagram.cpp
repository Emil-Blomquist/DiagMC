#include "FeynmanDiagram.h"

FeynmanDiagram::FeynmanDiagram (
  Vector3d externalMomentum,
  double length,
  double couplingConstant,
  double chemicalPotential
) {
  this->externalMomentum = externalMomentum;
  this->length = length;
  this->couplingConstant = couplingConstant;
  this->chemicalPotential = chemicalPotential;

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

double FeynmanDiagram::operator() () {
  double val;

  //
  // need to change when diagram length becomes variable
  //
  val = exp(this->chemicalPotential*this->length);

  for (auto g : this->Gs) {
    val *= (*g)();
  }

  for (auto d : this->Ds) {
    val *= (*d)(this->couplingConstant);
  }

  return val;
}

void FeynmanDiagram::save () {
  this->savedGs = this->Gs;
  this->savedDs = this->Ds;

  this->start->save();
  for (auto g : this->Gs) {
    g->save();
    g->end->save();
  }

  for (auto d : this->Ds) {
    d->save();
  }
}

void FeynmanDiagram::revert () {
  this->Gs = this->savedGs;
  this->Ds = this->savedDs;

  this->start->revert();
  for (auto g : this->Gs) {
    g->revert();
    g->end->revert();
  }

  for (auto d : this->Ds) {
    d->revert();
  }
}