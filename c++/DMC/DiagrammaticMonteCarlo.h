#ifndef DIAGRAMMATIC_MONTE_CARLO_H
#define DIAGRAMMATIC_MONTE_CARLO_H

#include <iostream>
#include <memory>
#include <vector>
#include <random>
#include <chrono> // seed mt
#include <math.h>
#include <boost/math/special_functions/erf.hpp>

#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

#include "../FeynmanDiagram/FeynmanDiagram.h"
#include "../Display/Display.h"

class DiagrammaticMonteCarlo {
    mt19937_64 mt;
    bool debug;
    int maxOrder;

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
};

#endif