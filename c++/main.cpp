#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <iomanip> // setprecision
#include <sstream> // stringstream

#include "DMC/DiagrammaticMonteCarlo.h"

using namespace std;

int main () {

  // parameters
  double
    p = 2,
    alpha = 2,
    mu = -2.2;

  // num iterations for each data point
  int N = 10000000;

  // external momentum
  Vector3d externalMomentum(p, 0, 0);

  // create vector with time data points
  vector<double> times(250, 0);
  for(int i = 0; i != times.size(); i++) {
    times[i] = 0.01 + i*5.0/(times.size() - 1);
  }

  // create date string to be written to file
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];
  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer, 80, "%Y-%m-%d %I:%M:%S", timeinfo);
  string dateString(buffer);

  // create file name
  stringstream stream;
  stream << fixed << setprecision(3)
         << "p="<< p << " a=" << alpha << " mu=" << mu;
  string fileName = stream.str();

  // write to file
  ofstream myfile;
  myfile.open("../data/" + fileName + ".txt", ios_base::app);
  myfile << "-------- " + dateString + " --------" << endl;
  myfile.close();

  // for each of our time data points
  for(int i = 0; i != times.size(); i++) {
    // initiate DMC
    DiagrammaticMonteCarlo DMC(externalMomentum, times[i], alpha, mu);

    // obtain g0 before we update the diagram
    double g0 = DMC.FD();

    // run the algorithm
    auto orderWeights = DMC.run(N, 7);

    // data to write to file
    myfile.open("../data/" + fileName + ".txt", ios_base::app);
    myfile << fixed << setprecision(7) << times[i] << " ";
    for (auto order : orderWeights) {
      myfile << fixed << setprecision(7) << g0*order/orderWeights[0] <<  " ";
    }
    myfile << fixed << setprecision(7) << g0/orderWeights[0] << endl;
    myfile.close();
  }
}