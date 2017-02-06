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


shared_ptr<Vertex> FeynmanDiagram::insertVertex (int gIndex, double dt) {
  if (dt <= 0) {
    cout << "ERROR at FeynmanDiagram::insertVertex: dt <= 0" << endl;
    return NULL;
  }

  shared_ptr<Electron> g = this->Gs[gIndex];

  double t = g->start->position + dt;
  double tmax = g->end->position;

  if (tmax <= t) {
    cout << "ERROR at FeynmanDiagram::insertVertex: tmax <= t" << endl;
    return NULL;
  }

  // create new vertex
  shared_ptr<Vertex> v(new Vertex(t));

  // create electronic propagator
  this->Gs.emplace(this->Gs.begin() + gIndex + 1, new Electron(externalMomentum));
  auto g1 = this->Gs[gIndex], g2 = this->Gs[gIndex + 1];

  // link them accordingly
  g2->setEnd(g1->end);
  g1->setEnd(v);
  g2->setStart(v);

  // return pointer to vertex
  return v;
}


shared_ptr<Phonon> FeynmanDiagram::addInternalPhonon (shared_ptr<Vertex> v1, shared_ptr<Vertex> v2, Vector3d Q, double theta, double phi) {
  // first make sure that v1 and v2 does not already are attached to a phonon
  if (v1->D[0] != NULL or v1->D[1] != NULL or v2->D[0] != NULL or v2->D[1] != NULL) {
    cout << "ERROR at FeynmanDiagram::addInternalPhonon: a phonon is already attached to either/both of these vertices" << endl;
    return NULL;
  }

  // create new phonon
  this->Ds.emplace_back(new Phonon(Q, theta, phi));

  // attach to vertices
  this->Ds.back()->setStart(v1);
  this->Ds.back()->setEnd(v2);

  // momentum conservation
  shared_ptr<Electron> g = v1->G[1];
  while (g->start != v2) {
    g->addMomentum(-Q);

    g = g->end->G[1];
  }

  return this->Ds.back();
}