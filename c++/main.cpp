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
  const double
    // maxMomenta = 3,
    // maxLength = 5*10,
    maxMomenta = 0,
    maxLength = 5,
    alpha = 2,
    mu = -2.2;

  const unsigned int
    // numIterations = 600*1000000,
    // numBins = 250*10,
    // maxOrder = 100000000;
    numIterations = 25000000,
    numBins = 250,
    maxOrder = 100000000;

  VectorXf momenta = VectorXf::LinSpaced(nprocs, 0, maxMomenta);
  Vector3d externalMomentum(momenta[myrank], 0, 0);

  double param = myrank;

  DiagrammaticMonteCarloV2 DMC(
    externalMomentum,
    maxLength,
    alpha,
    mu,
    numIterations,
    maxOrder,
    numBins,
    param
  );

  // do before MPI finalize
  printf("[Finished in %.2fs]\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

  MPI_Finalize();
  return 0;
}