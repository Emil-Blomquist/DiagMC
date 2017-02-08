#include "FeynmanDiagram.h"

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

void FeynmanDiagram::removeInternalPhonon (shared_ptr<Phonon> d) {
  // momentum conservation
  shared_ptr<Electron> g = d->start->G[1];
  while (g->start != d->end) {
    g->addMomentum(d->momentum);

    g = g->end->G[1];
  }

  // unlink phonon from vertices
  d->setStart(NULL);
  d->setEnd(NULL);

  //
  // IMPROVEMENT: we could implement a hash table for Ds, with the key being the adress to the Phonon instance.
  //

  // remove d from Ds
  for (auto i = this->Ds.begin(); i != this->Ds.end(); ++i) {
    if (*i == d) {
      this->Ds.erase(i);
      break;
    }
  }
}

void FeynmanDiagram::setInternalPhononMomentum (shared_ptr<Phonon> d, Vector3d Q) {
  Vector3d dQ = Q - d->momentum;

  // add momentum difference
  d->addMomentum(dQ);

  // momentum conservation
  shared_ptr<Electron> g = d->start->G[1];
  while (g->start != d->end) {
    g->addMomentum(-dQ);

    g = g->end->G[1];
  }
}

void FeynmanDiagram::setInternalPhononMomentumDirection (shared_ptr<Phonon> d, double theta, double phi) {
  d->setTheta(theta);
  d->setPhi(phi);
}