#ifndef FEYNMAN_DIAGRAM_H
#define FEYNMAN_DIAGRAM_H

#include <iostream>
#include <memory>
#include <vector>
#include <float.h>

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

#include "../Electron/Electron.h"
#include "../Phonon/Phonon.h"
#include "../Vertex/Vertex.h"

class FeynmanDiagram {
  private:
    unsigned int binaryElectronSearch(shared_ptr<Electron>, unsigned int lowerBound = 0);

  public:
    double length, couplingConstant, chemicalPotential;
    Vector3d externalMomentum;

    shared_ptr<Vertex> start, end;
    vector<shared_ptr<Electron> > Gs;
    vector<shared_ptr<Phonon> > Ds;

    FeynmanDiagram (Vector3d, double, double, double);

    shared_ptr<Vertex> insertVertex (int, double);
    shared_ptr<Phonon> addInternalPhonon (shared_ptr<Vertex>, shared_ptr<Vertex>, Vector3d, double, double, double);
    void
      removeVertex (shared_ptr<Vertex>),
      setVertexPosition (shared_ptr<Vertex>, double),
      swapPhonons (shared_ptr<Vertex>, shared_ptr<Vertex>),
      removeInternalPhonon (unsigned int),
      setInternalPhononMomentum (shared_ptr<Phonon>, Vector3d, double),
      setInternalPhononMomentumDirection (shared_ptr<Phonon>, double, double),
      setLength (double);

    double
      phononEnergy (double),
      operator() ();
};

#endif