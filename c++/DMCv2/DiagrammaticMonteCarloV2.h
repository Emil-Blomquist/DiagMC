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
#include <algorithm>

#include <fstream>
#include <string>
#include <ctime>
#include <iomanip> // setprecision
#include <sstream> // stringstream

#include <tuple> // return value

#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

#include "../FeynmanDiagram/FeynmanDiagram.h"
// #include "../Display/Display.h"

class DiagrammaticMonteCarloV2 {
  private:
    mt19937_64 mt;
    bool debug, loud;
    unsigned int maxOrder, numIterations, numBins;
    double maxLength, mu, param, lastKeyMargin;
    struct tm *timeinfo;

    vector<double> keys, bins;

    shared_ptr<Vertex> vertices2beRemoved [2];
    shared_ptr<Electron> electrons2beRemoved [2];
    shared_ptr<Phonon> phonon2beRemoved;

    double Udouble (double, double);
    int Uint (int, int);

    Vector3d calculateP0 (shared_ptr<Phonon>);
    Vector3d calculateQ (Vector3d, double, double, double);

    double shiftVertexPosition (double param = 1);
    double swapPhononConnections (double param = 1);
    double changeInternalPhononMomentumDirection (double param = 1);
    double changeInternalPhononMomentumMagnitude (double param = 1);
    double raiseOrder (double param = 1);
    double lowerOrder (double param = 1);
    double changeDiagramLength (double param = 1);

    void write2file (const unsigned int = 0);

  public:
    FeynmanDiagram FD;

    DiagrammaticMonteCarloV2 (
      Vector3d, 
      double,
      double,
      double,
      unsigned int,
      unsigned int,
      unsigned int,
      double param
    );
    
    void run ();
};

#endif