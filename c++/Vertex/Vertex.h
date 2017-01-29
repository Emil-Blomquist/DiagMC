#ifndef VERTEX_H
#define VERTEX_H

#include <iostream>

using namespace std;

class Propagator;

class Vertex {
  public:
    double position;

    Propagator * G [2];

    Vertex (double);

    void setIngoingG (Propagator *);
    void setOutgoingG (Propagator *);

    void print ();
};

#endif