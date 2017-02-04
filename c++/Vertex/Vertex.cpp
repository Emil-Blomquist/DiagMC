#include "Vertex.h"

Vertex::Vertex (double position) {
  this->position = position;

  this->G[0] = NULL;
  this->G[1] = NULL;
  this->D[0] = NULL;
  this->D[1] = NULL;
}

void Vertex::print () {
  cout << "Vertex: " << this << endl
       << "\tposition: " << position << endl
       << "\tG in: " << this->G[0] << endl
       << "\tG out: " << this->G[1] << endl
       << "\tD in: " << this->D[0] << endl
       << "\tD out: " << this->D[1] << endl;
}

// void Vertex::setPosition (double t) {
//   if (this->G[0] == NULL or this->G[1] == NULL) {
//     cout << "ERROR at Vertex::setPosition: cannot change position of first nor last vertex";
//     return;
//   }

//   double tmin = this->G[0]->start.position, tmax = this->G[1]->end.position;

// }