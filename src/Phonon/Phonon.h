#ifndef PHONON_H
#define PHONON_H

#include <iostream>
#include <memory>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

#include "../Vertex/Vertex.h"

class Phonon: public enable_shared_from_this<Phonon> {
  public:
    double q, theta, phi;

    Vector3d momentum;

    shared_ptr<Vertex> start, end;

    Phonon (Vector3d, double, double, double);

    void
      setStart (shared_ptr<Vertex>),
      setEnd (shared_ptr<Vertex>),
      setMomentum (Vector3d, double),
      setTheta (double),
      setPhi (double);

    double operator() (double, double);
};

#endif