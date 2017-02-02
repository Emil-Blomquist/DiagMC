#ifndef FEYNMAN_DIAGRAM_H
#define FEYNMAN_DIAGRAM_H

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

#include "../Propagator/Electron.h"
#include "../Propagator/Phonon.h"
#include "../Vertex/Vertex.h"

class FeynmanDiagram {
  public:
    double length, couplingConstant, chemicalPotential;
    Vector3d externalMomentum;

    vector<Electron*> Gs;
    vector<Phonon*> Ds;

    Vertex *start, *end;

    FeynmanDiagram (Vector3d, double, double, double);
};

#endif