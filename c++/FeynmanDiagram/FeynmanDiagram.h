#ifndef FEYNMAN_DIAGRAM_H
#define FEYNMAN_DIAGRAM_H

#include <iostream>
#include <memory>
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

    shared_ptr<Vertex> start, end;

    vector<shared_ptr<Electron>> Gs;
    vector<shared_ptr<Phonon>> Ds;

    FeynmanDiagram (Vector3d, double, double, double);

    void print ();

    shared_ptr<Vertex> insertVertex (int, double);
    shared_ptr<Phonon> addInternalPhonon (shared_ptr<Vertex>, shared_ptr<Vertex>, Vector3d, double, double);
};

#endif