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
    unsigned int numBins, n00;
    unsigned long int numIterations;
    double maxLength, mu, alpha, param, lastKeyMargin;
    struct tm *timeinfo;
    char **argv;

    vector<double> keys;
    vector<int> bins;

    shared_ptr<Vertex> vertices2beRemoved [2];
    shared_ptr<Electron> electrons2beRemoved [2];
    shared_ptr<Phonon> phonon2beRemoved;

    double Udouble (double, double);
    int Uint (int, int);

    Vector3d
      calculateMeanP (shared_ptr<Vertex>, shared_ptr<Vertex>),
      calculateP0 (shared_ptr<Phonon>);

    void
      shiftVertexPosition (double param = 1),
      swapPhononConnections (double param = 1),
      changeInternalPhononMomentumDirection (double param = 1),
      changeInternalPhononMomentumMagnitude (double param = 1),
      raiseOrder (double param = 1),
      lowerOrder (double param = 1),
      changeDiagramLength (double param = 1);

    void write2file (const unsigned long int = 0);

  public:
    FeynmanDiagram FD;

    DiagrammaticMonteCarloV2 (
      Vector3d, 
      double,
      double,
      double,
      unsigned long int,
      unsigned int,
      double param,
      char **argv
    );
    
    void run ();
};

#endif