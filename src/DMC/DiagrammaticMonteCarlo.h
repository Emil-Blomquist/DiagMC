#ifndef DIAGRAMMATIC_MONTE_CARLO_H
#define DIAGRAMMATIC_MONTE_CARLO_H

#include <iostream>
#include <memory>
#include <vector>
#include <random>
#include <chrono> // seed mt
#include <math.h>
#include <boost/math/special_functions/erf.hpp> // inverse error function
#include <limits>
#include <algorithm>
#include <map>

#include <fstream>
#include <string>
#include <ctime>
#include <iomanip> // setprecision
#include <sstream> // stringstream

#include <tuple> // return value

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

#include "../FeynmanDiagram/FeynmanDiagram.h"
#include "../Display/Display.h"

class DiagrammaticMonteCarlo {
  private:
    mt19937_64 mt;
    bool debug, loud, externalLegs, irreducibleDiagrams, fixedExternalMomentum;
    unsigned long int numIterations, N0;
    double maxLength, mu, alpha, param, dt, dp, maxMomenta;
    struct tm *timeinfo;
    char **argv;

    vector<unsigned long int> bins, bins0;

    Array<unsigned long int, Dynamic, Dynamic> hist;

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
      changeDiagramLength (double param = 1),
      changeDiagramLengthComplex (double param = 1),
      changeExternalMomentumMagnitude (double param = 1);

    void
      write2file (const unsigned long int = 0),
      checkAcceptanceRatio (double, string);

  public:
    FeynmanDiagram FD;

    DiagrammaticMonteCarlo (
      Vector3d, 
      double,
      double,
      double,
      unsigned long int,
      // unsigned int,
      double param,
      char **argv
    );
    
    void run ();
};

#endif