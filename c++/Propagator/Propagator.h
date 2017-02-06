#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include <iostream>
#include <memory>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

#include "../Vertex/Vertex.h"

class Propagator {
  public:
    int type;
    bool removed;

    Vector3d momentum;

    shared_ptr<Vertex> start, end;

    Propagator (int, Vector3d);

    void setMomentum (Vector3d);
    void addMomentum (Vector3d);
};

#endif