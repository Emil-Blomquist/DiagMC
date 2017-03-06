#include <iostream>
#include <time.h> // execution time
#include <mpi.h> // MPI

#include "DMCv2/DiagrammaticMonteCarloV2.h"

using namespace std;

int main () {
  // MPI stuff
  int myrank, nprocs;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // do after MPI stuff
  clock_t tStart = clock();

  // parameters
  // const double
  //   minMomenta = 1.4,
  //   maxMomenta = 1.8,
  //   maxLength = 30,
  //   // maxMomenta = 0,
  //   // maxLength = 5*2,
  //   alpha = 1,
  //   mu = -0.3;

  // const unsigned int
  //   numIterations = 3000000000,
  //   numBins = 250*6;
  //   // numIterations = 10000000*2,
  //   // numBins = 250*2;


  const double
    minMomenta = 0,
    maxMomenta = 0,
    maxLength = 10,
    alpha = 2,
    mu = -2.2;

  const unsigned int
    numIterations = 10000000*2,
    numBins = 250*2;

  VectorXf momenta = VectorXf::LinSpaced(nprocs, minMomenta, maxMomenta);
  Vector3d externalMomentum(momenta[myrank], 0, 0);

  double param = myrank;

  DiagrammaticMonteCarloV2 DMC(
    externalMomentum,
    maxLength,
    alpha,
    mu,
    numIterations,
    numBins,
    param
  );


  // do before MPI finalize
  printf("[Finished in %.2fs]\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

  MPI_Finalize();
  return 0;
}