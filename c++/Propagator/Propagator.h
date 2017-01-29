#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

class Vertex;

class Propagator {
    Vector3d savedMomentum;

  public:
    int type;

    Vector3d momentum;
    double theta, phi;

    Vertex * start, * end;

    Propagator (int, Vector3d, double, double);
    // ~Propagator();

    void setMomentum (Vector3d);
    void setStart (Vertex *);
    void setEnd (Vertex *);
    
    void print ();
};

#endif