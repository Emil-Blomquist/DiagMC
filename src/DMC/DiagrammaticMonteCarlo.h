#ifndef DIAGRAMMATIC_MONTE_CARLO_H
#define DIAGRAMMATIC_MONTE_CARLO_H

#include <mpi.h> // MPI

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
#include <unsupported/Eigen/FFT>

using namespace Eigen;
using namespace std;

#include "../FeynmanDiagram/FeynmanDiagram.h"
// #include "../Display/Display.h"

class DiagrammaticMonteCarlo {
  private:
    mt19937_64 mt;
    bool debug, loud, externalLegs, reducibleDiagrams, fixedExternalMomentum, Dyson, skeletonDiagrams, bold;
    unsigned long long int numIterations, N0, numMCIterations, untilStart, currItr;
    double maxLength, mu, alpha, param, dt, dp, maxMomenta;
    struct tm *timeinfo;
    char **argv;
    unsigned int Np, Nt, minDiagramOrder, maxDiagramOrder, numBoldIterations, boldIteration, MCvsDMCboundary;

    unsigned int savePeriod;

    int worldRank, worldSize;

    Array<unsigned long long int, Dynamic, Dynamic> hist;
    Array<double, Dynamic, Dynamic> S, dG, dE, S1mc;

    shared_ptr<Vertex> vertices2beRemoved [2];
    shared_ptr<Electron> electrons2beRemoved [2];
    shared_ptr<Phonon> phonon2beRemoved;

    Vector3d initialExternalMomentum;

    // structure to be sent to rank 0 when calculating S1
    struct KeyValue {
      unsigned int pi;
      unsigned int ti;
      double val;
    };

    double
      dEOf (shared_ptr<Electron>),
      dEOf (double, double),
      additionalPhase (shared_ptr<Electron>),
      additionalPhase (Vector3d, double),
      additionalPhase (double, double),
      evaluateDiagram (),
      Udouble (double, double),
      Ndouble (double);

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
      BOLDchangeDiagramLength (double param = 1),
      BOLDchangeDiagramLengthComplex (double param = 1),
      BOLDraiseOrder (double param = 1),
      BOLDlowerOrder (double param = 1),
      BOLDshiftVertexPosition (double param = 1);

    void
      write2file (
        Array<unsigned long long int, Dynamic, Dynamic>&,
        unsigned long long int&,
        unsigned long long int&),
      doDyson (ArrayXXd&, ArrayXXd&),
      normalizeHistogram (
        Array<unsigned long long int, Dynamic, Dynamic>&,
        unsigned long long int&,
        ArrayXXd&),
      calculateEnergyDiff (ArrayXXd&, ArrayXXd&),
      checkAcceptanceRatio (double, string);

    void
      firstOrderSelfEnergyMC (),
      appendKeyValue (vector<KeyValue>&, unsigned int), // <- gather MC calculation
      sumHistograms (
        Array<unsigned long long int, Dynamic, Dynamic>&,
        unsigned long long int&,
        unsigned long long int&),
      importG (string);

    double firstOrderSelfEnergyMC (double, Vector3d);

  public:
    FeynmanDiagram FD;

    DiagrammaticMonteCarlo (
      Vector3d, 
      double,
      double,
      double,
      unsigned long long int,
      // unsigned int,
      double param,
      char **argv
    );
    
    void run ();
};

#endif