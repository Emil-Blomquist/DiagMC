#ifndef ELECTRON_H
#define ELECTRON_H

#include <iostream>
#include <memory>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

#include "../Vertex/Vertex.h"

class Electron: public enable_shared_from_this<Electron> {
  public:
    double p;
    Vector3d momentum;

    shared_ptr<Vertex> start, end;

    Electron (Vector3d);

    void
      setStart (shared_ptr<Vertex>),
      setEnd (shared_ptr<Vertex>),
      addMomentum (Vector3d);
    
    double operator() (double, double dE = 0);

    static double value (double, double, double, double dE = 0);
};

#endif