#include "FeynmanDiagram.h"

shared_ptr<Phonon> FeynmanDiagram::addInternalPhonon (
  shared_ptr<Vertex> v1,
  shared_ptr<Vertex> v2,
  Vector3d Q, 
  double q,
  double theta,
  double phi
) {
  // first make sure that v1 and v2 are not already attached to a phonon
  if (v1->D[0] != NULL or v1->D[1] != NULL or v2->D[0] != NULL or v2->D[1] != NULL) {
    cout << "ERROR at FeynmanDiagram::addInternalPhonon: a phonon is already attached to either/both of these vertices" << endl;
    return NULL;
  }

  // create new phonon
  this->Ds.emplace_back(new Phonon(Q, q, theta, phi));

  // attach to vertices
  this->Ds.back()->setStart(v1);
  this->Ds.back()->setEnd(v2);

  // momentum conservation
  shared_ptr<Electron> g = v1->G[1];
  do {
    // remove from hash table
    this->removeFromHashTable(g);  

    g->addMomentum(-Q);

    // reinsert into hash table
    this->insertIntoHashTable(g);
  } while (g->end != v2 && (g = g->end->G[1]));

  return this->Ds.back();
}

void FeynmanDiagram::removeInternalPhonon (unsigned int phononIndex) {
  shared_ptr<Phonon> d = this->Ds[phononIndex];

  // momentum conservation
  shared_ptr<Electron> g = d->start->G[1];

  do {
    // remove from hash table
    this->removeFromHashTable(g);  

    g->addMomentum(d->momentum);

    // reinsert into hash table
    this->insertIntoHashTable(g);
  } while (g->end != d->end && (g = g->end->G[1]));

  // unlink phonon from vertices
  d->setStart(NULL);
  d->setEnd(NULL);

  // remove d from Ds
  this->Ds.erase(this->Ds.begin() + phononIndex);
}

void FeynmanDiagram::setInternalPhononMomentum (shared_ptr<Phonon> d, Vector3d Q, double q) {
  // obviously need to be calculated before adding the momentum
  Vector3d dQ = Q - d->momentum;

  // add momentum difference
  d->setMomentum(Q, q);

  // momentum conservation
  shared_ptr<Electron> g = d->start->G[1];
  do {
    // remove from hash table
    this->removeFromHashTable(g);  

    g->addMomentum(-dQ);
    
    // reinsert into hash table
    this->insertIntoHashTable(g);
  } while (g->end != d->end && (g = g->end->G[1]));
}

void FeynmanDiagram::setInternalPhononMomentumDirection (shared_ptr<Phonon> d, double theta, double phi) {
  d->setTheta(theta);
  d->setPhi(phi);
}