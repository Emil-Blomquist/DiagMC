#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

// needed since we will call methods
#include "../Vertex/Vertex.h"

class Propagator {
  public:
    int type;

    Vector3d momentum;

    Vertex *start, *end;

    Propagator (int, Vector3d);
    // ~Propagator();

    void setMomentum (Vector3d);
    void setStart (Vertex *);
    void setEnd (Vertex *);
    
    void print ();
};

#endif