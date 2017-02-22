#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <iomanip> // setprecision
#include <sstream> // stringstream
#include <thread>

#include "DMC/DiagrammaticMonteCarlo.h"

using namespace std;


void runParallel (
  double& state,
  const VectorXf& times,
  const double p,
  const double alpha,
  const double mu,
  const int numIterations,
  const int maxOrder
) {
  // external momentum
  Vector3d externalMomentum(p, 0, 0);

  // create date and time string
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];
  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer, 80, "%Y-%m-%d %I:%M:%S", timeinfo);
  string dateAndTimeString(buffer);

  // create file name
  stringstream stream;
  stream << fixed << setprecision(3)
         << "p="<< p << " a=" << alpha << " mu=" << mu << " N=" << numIterations << " time=" << dateAndTimeString;
  string fileName = stream.str();

  // write to file
  ofstream myfile;
  myfile.open("../data/" + fileName + ".txt", ios_base::app);
  myfile << "-------- " + fileName + " --------" << endl;
  myfile.close();

  // for each of our time data points
  for(int i = 0; i != times.size(); i++) {
    // initiate DMC
    DiagrammaticMonteCarlo DMC(externalMomentum, times[i], alpha, mu);

    // obtain g0 before we update the diagram
    double g0 = DMC.FD();

    // run the algorithm
    auto orderWeights = DMC.run(numIterations, maxOrder);

    // data to write to file
    myfile.open("../data/" + fileName + ".txt", ios_base::app);
    myfile << fixed << setprecision(7) << times[i] << " ";
    for (auto order : orderWeights) {
      myfile << fixed << setprecision(7) << g0*order/orderWeights[0] <<  " ";
    }
    myfile << fixed << setprecision(7) << g0/orderWeights[0] << endl;
    myfile.close();

    // increase state value
    state = (i + 1.0)/times.size();
  }

  // in case we have a rounding error
  state = 1;
}








int main () {
  // parameters
  const double
    pMax = 5,
    tMax = 5,
    alpha = 2,
    mu = -2.2;

  // num iterations for each data point
  const int numIterations = 10000000;

  // maximum diagram order considered
  const int maxOrder = 7;

  // create vector with time data points
  VectorXf times = VectorXf::LinSpaced(250, 0.01, tMax + 0.01);

  // number of threads able to run parallel
  unsigned int numThreadsToUse = thread::hardware_concurrency()/2;

  // create vector with time data points
  VectorXf momenta = VectorXf::LinSpaced(numThreadsToUse, 0, pMax);

  // this will store the state of each parallel running program
  vector<double> state(numThreadsToUse, 0);

  // run parallel code
  vector<thread> threads;
  for (int i = 0; i < numThreadsToUse; ++i) {
    threads.push_back(thread(runParallel, ref(state[i]), ref(times), momenta[i], alpha, mu, numIterations, maxOrder));
  }

  // display state of simulations
  while (true) {
    cout << fixed << setprecision(3);
    for(auto &c : state){
      cout << c << "\t";
    }
    cout << endl;

    // determin whether we are done
    bool wereDone = true;    
    for(auto &c : state){
      if (c < 1) {
        wereDone = false;
        break;
      }
    }

    // if we are done
    if (wereDone) {
      break;
    }

    // sleep for 15s
    usleep(1000000*15);
  }

  // awaiting parallel code
  for(auto &t : threads){
    t.join();
  }

  cout << "CALCULATIONS COMPLETE" << endl;
  return 0;
}