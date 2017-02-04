#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

#include "../Vertex/Vertex.h"

class Propagator {
  public:
    int type;
    bool dummy;

    Vector3d momentum;

    Vertex *start, *end;

    Propagator ();
    Propagator (int, Vector3d);

    void setMomentum (Vector3d);
    void addMomentum (Vector3d);
};

#endif