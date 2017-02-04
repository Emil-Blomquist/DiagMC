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

  // create first and last vertex
  this->Vs.emplace_back(0);
  this->Vs.emplace_back(length);

  // create a electronic propagator
  this->Gs.emplace_back(externalMomentum);

  // link propagator to vertices
  this->Gs.begin()->setStart(&this->Vs.front());
  this->Gs.begin()->setEnd(&this->Vs.back());
}

void FeynmanDiagram::print () {
  // print vertexes
  for (list<Vertex>::iterator i = this->Vs.begin(); i != this->Vs.end(); ++i) {
    i->print();
  }

  // print electronic propagators
  for (list<Electron>::iterator j = this->Gs.begin(); j != this->Gs.end(); ++j) {
    j->print();
  }

  // print phonon propagators
  for (list<Phonon>::iterator k = this->Ds.begin(); k != this->Ds.end(); ++k) {
    k->print();
  }
}


Vertex *FeynmanDiagram::insertVertex (list<Vertex>::iterator i, double dt) {
  if (dt <= 0) {
    cout << "ERROR at FeynmanDiagram::insertVertex: dt <= 0" << endl;
    return NULL;
  }

  double t = i->position + dt;
  double tmax = (++i)->position;

  if (tmax <= t) {
    cout << "ERROR at FeynmanDiagram::insertVertex: tmax <= t" << endl;
    return NULL;
  }

  // create new vertex
  Vertex *v;
  this->Vs.emplace(i, t);
  v = &(*--i);

  // create new electronic propagator
  Electron *g1, *g2;
  g1 = (--i)->G[1];
  this->Gs.emplace_back(g1->momentum);
  g2 = &this->Gs.back();

  // link them accordingly
  g2->setEnd(g1->end);
  g1->setEnd(v);
  g2->setStart(v);

  // return pointer to vertex
  return v;
}


Phonon *FeynmanDiagram::addInternalPhonon (list<Vertex>::iterator v1, list<Vertex>::iterator v2, Vector3d Q, double theta, double phi) {
  // first make sure that v1 and v2 does not already are attached to a phonon
  if (v1->D[0] != NULL or v1->D[1] != NULL or v2->D[0] != NULL or v2->D[1] != NULL) {
    cout << "ERROR at FeynmanDiagram::addInternalPhonon: a phonon is already attached to either/both of these vertices" << endl;
    return NULL;
  }

  // create new phonon
  Phonon *d;
  this->Ds.emplace_back(Q, theta, phi);
  d = &this->Ds.back();

  // attach to vertices
  d->setStart(&(*v1));
  d->setEnd(&(*v2));

  //
  // conserve total momentum
  //
  for (list<Vertex>::iterator i = v1; i != v2; ++i) {
    i->G[1]->addMomentum(-Q);
  }

  return d;
}