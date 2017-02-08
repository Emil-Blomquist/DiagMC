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

    // vertex methods
    shared_ptr<Vertex> insertVertex (int, double);
    void removeVertex (shared_ptr<Vertex>);
    void setVertexPosition (shared_ptr<Vertex>, double);
    void swapPhonons (shared_ptr<Vertex>, shared_ptr<Vertex>);
    // phonon methods
    shared_ptr<Phonon> addInternalPhonon (shared_ptr<Vertex>, shared_ptr<Vertex>, Vector3d, double, double);
    void removeInternalPhonon (shared_ptr<Phonon>);
    void setInternalPhononMomentum (shared_ptr<Phonon>, Vector3d);
    void setInternalPhononMomentumDirection (shared_ptr<Phonon>, double, double);

};

#endif