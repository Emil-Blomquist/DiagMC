#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <iostream>
#include <memory>
#include <vector>
#include <random>
#include <chrono> // seed mt
#include <math.h>

#include <fstream>
#include <string>
#include <ctime>
#include <iomanip> // setprecision
#include <sstream> // stringstream

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

class MonteCarlo {
  private:
    mt19937_64 mt;
    struct tm *timeinfo;
    string path;

    char **argv;

    bool externalLegs, irreducibleDiagrams;

    Vector3d externalMomentum;
    unsigned long int numIterations;
    double mu, alpha, tStart, tEnd, dtOut, dt, dp, maxLength, maxMomenta;
    vector<double> times, values;

    unsigned int maxOrder;

    Array<double, Dynamic, Dynamic> Gfull, dG, dE;

    double
      Udouble (double, double),
      Ndouble (double);

    void
      write2file (),
      run (),
      diagramOrder1 (double, unsigned int),
      diagramOrder2 (double, unsigned int),
      importG (string),
      calculateEnergyDiff ();

    double
      phononDispersionRelation (Vector3d),
      G (Vector3d, double, double),
      D (Vector3d, double, double, double);

  public:
    MonteCarlo (
      Vector3d,
      double,
      double,
      unsigned long int,
      double,
      double,
      char**
    );
};

#endif