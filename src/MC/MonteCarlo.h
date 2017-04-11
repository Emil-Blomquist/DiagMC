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

    bool externalLegs, irreducibleDiagrams;

    Vector3d externalMomentum;
    unsigned long int numIterations;
    double mu, alpha, tMax;
    vector<double> times, values;

    double
      Udouble (double, double),
      Ndouble (double);

    void
      write2file (),
      run (),
      diagramOrder1 (double, unsigned int),
      diagramOrder2 (double, unsigned int);

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
      unsigned int,
      char **argv
    );
};

#endif