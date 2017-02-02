#include "Vertex.h"

Vertex::Vertex (double _position) {
  position = _position;

  G[0] = NULL; G[1] = NULL;
  D[0] = NULL; D[1] = NULL;
}

void Vertex::print () {
  cout << "Vertex:" << endl
       << "\tposition: " << position << endl
       << "\tingoing g: " << G[0] << endl
       << "\toutgoing g: " << G[1] << endl
       << "\tingoing d:" << D[0] << endl
       << "\toutgoing d:" << D[1] << endl;

}