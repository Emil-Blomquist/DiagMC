#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include <iostream>
#include <memory>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

#include "../Vertex/Vertex.h"

class Propagator {
    Vector3d savedMomentum;
    shared_ptr<Vertex> savedStart, savedEnd;

  public:
    int type;

    Vector3d momentum;

    shared_ptr<Vertex> start, end;

    Propagator (int, Vector3d);

    void setMomentum (Vector3d);
    void addMomentum (Vector3d);

    void save ();
    void revert ();

    void unlink ();
};

#endif