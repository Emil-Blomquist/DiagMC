#include "Vertex.h"
#include "../Propagator/Propagator.h"

Vertex::Vertex (double _position) {
  position = _position;

  G[0] = NULL; G[1] = NULL;
}


void Vertex::setIngoingG (Propagator * g) {
  G[0] = g;
}

void Vertex::setOutgoingG (Propagator * g) {
  G[1] = g;
}

void Vertex::print () {
  cout << "Vertex:" << endl
       << "\tposition: " << position << endl
       << "\tingoing g: " << G[0] << endl
       << "\toutgoing g: " << G[1] << endl;

}