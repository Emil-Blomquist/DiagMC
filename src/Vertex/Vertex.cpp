#include "Vertex.h"

Vertex::Vertex (double position) {
  this->position = position;
}


void Vertex::setG (int i, shared_ptr<Electron> g) {
  if (g == NULL) {
    this->G[i].reset();
  } else {
    this->G[i] = g;
  }
}

void Vertex::setD (int i, shared_ptr<Phonon> d) {
  if (d == NULL) {
    this->D[i].reset();
  } else {
    this->D[i] = d;
  }
}

void Vertex::setPosition (double t) {
  if ( ! isfinite(t)) {
    cout << "Phonon::setPosition t=" << t << endl;
  }

  this->position = t;
}