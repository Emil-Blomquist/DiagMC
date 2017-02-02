#ifndef VERTEX_H
#define VERTEX_H

#include <iostream>

using namespace std;

class Propagator;

class Vertex {
  public:
    double position;

    Propagator *G[2];
    Propagator *D[2];

    Vertex (double);

    void setIngoingG (Propagator *g) { G[0] = g; }
    void setOutgoingG (Propagator *g) { G[1] = g; }
    void setIngoingD (Propagator *d) { D[0] = d; }
    void setOutgoingD (Propagator *d) { D[1] = d; }


    void print ();
};

#endif