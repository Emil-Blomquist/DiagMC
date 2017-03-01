#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <iomanip> // setprecision
#include <sstream> // stringstream
#include <thread>

#include <time.h> // execution time

#include <mpi.h> // MPI

#include "DMCv2/DiagrammaticMonteCarloV2.h"

// using namespace std;


void runParallel (
  const double p,
  const double maxLength,
  const double alpha,
  const double mu,
  const int numIterations,
  const int numBins,
  const int maxOrder,
  const int processId,
  const int numProcesses
) {
  // external momentum
  Vector3d externalMomentum(p, 0, 0);

  // create date and time string
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];
  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", timeinfo);
  string dateAndTimeString(buffer);

  VectorXf params = VectorXf::LinSpaced(numProcesses, 0.1, 4);
  double param = params[processId];

  // create file name
  stringstream stream;
  stream << fixed << setprecision(3)
         << "p="<< p
         << " tmax=" << maxLength
         << " a=" << alpha
         << " mu=" << mu
         << " N=" << numIterations
         << " date=" << dateAndTimeString
         << " id=" << processId
         << " param=" << param;
  string fileName = stream.str();

  // write to file
  ofstream myfile;
  myfile.open("../data/" + fileName + ".txt", ios_base::app);
  myfile << "-------- " + fileName + " --------" << endl;
  myfile.close();

  // for each of our time data points
  DiagrammaticMonteCarloV2 DMC(externalMomentum, maxLength, alpha, mu, param);

  // run the algorithm
  auto result = DMC.run(numIterations, maxOrder, numBins);

  auto
    keys = get<0>(result),
    bins = get<1>(result);

  // data to write to file
  myfile.open("../data/" + fileName + ".txt", ios_base::app);
  for (auto key : keys) {
    myfile << fixed << setprecision(7) << key <<  " ";
  }
  myfile << "\n";
  for (auto bin : bins) {
    myfile << fixed << setprecision(7) << bin <<  " ";
  }
  myfile.close();
}



int main () {
  // MPI stuff
  int myrank, nprocs;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // do after MPI stuff
  clock_t tStart = clock();

  // parameters
  const double
    momentum = 0,
    maxLength = 5,
    alpha = 2,
    mu = -2.2;

  const unsigned int
    numIterations = 25000000,
    numBins = 250,
    maxOrder = 7;

  runParallel(momentum, maxLength, alpha, mu, numIterations, numBins, maxOrder, myrank, nprocs);

  // do before MPI finalize
  printf("[Finished in %.2fs]\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

  MPI_Finalize();
  return 0;
}