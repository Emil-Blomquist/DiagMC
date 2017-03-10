#include <iostream>
#include <time.h> // execution time
// #include <mpi.h> // MPI

#include "DMCv2/DiagrammaticMonteCarloV2.h"

using namespace std;

int main () {
  clock_t tStart = clock();

  // parameters
  // const double
  //   minMomenta = 1.8,
  //   maxMomenta = 1.8,
  //   maxLength = 30,
  //   alpha = 1,
  //   mu = 0.03;

  // const unsigned int
    // numIterations = 3000000000,
  //   numBins = 250*6;


  const double
    minMomenta = 1.8,
    maxMomenta = 1.8,
    maxLength = 10,
    alpha = 1,
    mu = 0.03;

  const unsigned int
    // numIterations = 50000000*2,
    numIterations = 4000000000,
    numBins = 250*2;

  Vector3d externalMomentum(0, 0, minMomenta);

  double param = 0;

  // DiagrammaticMonteCarloV2 DMC(
  //   externalMomentum,
  //   maxLength,
  //   alpha,
  //   mu,
  //   numIterations,
  //   numBins,
  //   param
  // );



  printf("[Finished in %.2fs]\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  return 0;
}



// int myrank, nprocs;
// MPI_Init(NULL, NULL);
// MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
// MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

// ... code ...

// MPI_Finalize();