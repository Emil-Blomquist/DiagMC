#include "FeynmanDiagram.h"

shared_ptr<Vertex> FeynmanDiagram::insertVertex (int gIndex, double dt) {
  if ( ! isfinite(dt)) {
    cout << "FeynmanDiagram::insertVertex dt=" << dt << endl;
  }

  if (dt <= 0) {
    // cout << "ERROR at FeynmanDiagram::insertVertex: dt <= 0" << endl;
    // return NULL;
    dt = DBL_EPSILON;
  }

  shared_ptr<Electron> g = this->Gs[gIndex];

  double t = g->start->position + dt;
  double tmax = g->end->position;

  if (tmax <= t) {
    // cout << "ERROR at FeynmanDiagram::insertVertex: tmax <= t" << endl;
    // return NULL;
    t = tmax - DBL_EPSILON;
  }

  // create new vertex
  shared_ptr<Vertex> v(new Vertex(t));

  // create electronic propagator
  this->Gs.emplace(this->Gs.begin() + gIndex + 1, new Electron(g->momentum));
  auto g1 = this->Gs[gIndex], g2 = this->Gs[gIndex + 1];

  // link them accordingly
  g2->setEnd(g1->end);
  g1->setEnd(v);
  g2->setStart(v);

  // return pointer to vertex
  return v;
}

void FeynmanDiagram::removeVertex (shared_ptr<Vertex> v) {
  if (v->D[0] || v->D[1]) {
    cout << "ERROR at FeynmanDiagram::removeVertex: Phonon line still attatched to vertex" << endl;
    return;
  }

  if (v == this->start || v == this->end) {
    cout << "ERROR at FeynmanDiagram::removeVertex: Cannot remove the first nor the last vertex" << endl;
    return;
  }

  // temporary save
  shared_ptr<Electron> g = v->G[1];

  // obtain index before we unlink
  unsigned int i = this->binaryElectronSearch(g);

  // relink
  v->G[0]->setEnd(g->end);

  // unlink
  g->setStart(NULL);

  // remove g from Gs
  this->Gs.erase(this->Gs.begin() + i);   
}

void FeynmanDiagram::setVertexPosition (shared_ptr<Vertex> v, double t) {
  if (v == this->start || v == this->end) {
    cout << "ERROR at FeynmanDiagram::setVertexPosition: Not allowed to change position of the first nor the last vertex" << endl;
    return;
  }

  if ( ! isfinite(t)) {
    cout << "FeynmanDiagram::setVertexPosition t=" << t << endl;
  }

  // if ( ! (v->G[0]->start->position < t && t < v->G[1]->end->position)) {
  //   cout << "ERROR at FeynmanDiagram::setVertexPosition: New vertex time invalid" << endl;
  //   return;
  // }

  if (v->G[0]->start->position >= t) {
    t = v->G[0]->start->position + DBL_EPSILON;
  } else if (t >= v->G[1]->end->position) {
    t = v->G[1]->end->position - DBL_EPSILON;
  }

  v->setPosition(t);
}

void FeynmanDiagram::swapPhonons (shared_ptr<Vertex> v1, shared_ptr<Vertex> v2) {
  // cannot be first or last vertex
  if (v1 == this->start || v2 == this->start || v1 == this->end || v2 == this->end) {
    cout << "ERROR at FeynmanDiagram::swapPhonons: Not allowed to swap phonon with the first nor the last vertex" << endl;
    return;
  }

  // must be neighbouring vertices and sent here in proper order
  if (v1->G[1]->end != v2) {
    cout << "ERROR at FeynmanDiagram::swapPhonons: Vertices not correctly ordered or not neighbouring ones" << endl;
    return;
  }

  // we do not allow for vertices to be connected to the same phonon
  if (v1->D[1] && v1->D[1]->end == v2) {
    cout << "ERROR at FeynmanDiagram::swapPhonons: The vertices must not connect to the same phonon" << endl;
    return;
  }

  shared_ptr<Phonon> d1, d2;
  d1 = (v1->D[0]) ? v1->D[0] : v1->D[1];
  d2 = (v2->D[0]) ? v2->D[0] : v2->D[1];

  Vector3d dP1, dP2, dP;
  dP1 = d1->momentum * ((v1->D[0]) ? -1 : 1);
  dP2 = d2->momentum * ((v2->D[1]) ? -1 : 1);
  dP = dP1 + dP2;

  // momentum conservation
  v1->G[1]->addMomentum(dP);

  // vector of pointers to member function of Phonon
  vector<void (Phonon::*)(shared_ptr<Vertex>)> startOrEnd = {&Phonon::setStart, &Phonon::setEnd};

  // set start or set end
  int i1, i2;
  i1 = (d1->end == v1) ? 1 : 0;
  i2 = (d2->end == v2) ? 1 : 0;

  // relink
  (*d1.*startOrEnd[i1])(v2);
  (*d2.*startOrEnd[i2])(v1);
}