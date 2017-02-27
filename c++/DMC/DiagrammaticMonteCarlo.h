#ifndef DIAGRAMMATIC_MONTE_CARLO_H
#define DIAGRAMMATIC_MONTE_CARLO_H

#include <iostream>
#include <memory>
#include <vector>
#include <random>
#include <chrono> // seed mt
#include <math.h>
#include <boost/math/special_functions/erf.hpp>
#include <limits>

#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

#include "../FeynmanDiagram/FeynmanDiagram.h"
// #include "../Display/Display.h"

class DiagrammaticMonteCarlo {
    mt19937_64 mt;
    bool debug, loud;
    int maxOrder, numIterations;

    shared_ptr<Vertex> vertices2beRemoved [2];
    shared_ptr<Electron> electrons2beRemoved [2];
    shared_ptr<Phonon> phonon2beRemoved;

    double Udouble (double, double);
    int Uint (int, int);

    Vector3d calculateP0 (shared_ptr<Phonon>);
    Vector3d calculateQ (Vector3d, double, double, double);

    double shiftVertexPosition ();
    double swapPhononConnections ();
    double changeInternalPhononMomentumDirection ();
    double changeInternalPhononMomentumMagnitude ();
    double raiseOrder ();
    double lowerOrder ();

  public:
    FeynmanDiagram FD;

    DiagrammaticMonteCarlo (Vector3d, double, double, double);
    
    vector<double> run (int, int);
};

#endif