#include <mpi.h> // MPI
#include <iostream>
#include <time.h> // execution time

#include "DMC/DiagrammaticMonteCarlo.h"
#include "MC/MonteCarlo.h"

using namespace std;

int main (int argc, char **argv) {
  clock_t tStart = clock();

  // default parameters
  double
    momenta = 0,
    alpha = 1,
    mu = -1.2,
    maxLength = 30,
    param = 0;

  unsigned long long int numIterations = 4000000000;

  // input parameters
  int i = 1;
  while (i < argc) {
    if (strcmp(argv[i], "-p") == 0) {
      momenta = stod(argv[i + 1]);
      i += 2;
    } else if (strcmp(argv[i], "-a") == 0) {
      alpha = stod(argv[i + 1]);
      i += 2;
    } else if (strcmp(argv[i], "-mu") == 0) {
      mu = stod(argv[i + 1]);
      i += 2;
    } else if (strcmp(argv[i], "-t") == 0) {
      maxLength = stod(argv[i + 1]);
      i += 2;
    } else if (strcmp(argv[i], "-N") == 0) {
      numIterations = stol(argv[i + 1]);
      i += 2;
    } else {
      i++;
    }
  }



  int worldRank, worldSize;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  if (worldRank == 0) {
    cout << endl
         << "1) implement numSecs instead of numItr" << endl
         << endl;
  }


  Vector3d externalMomentum(0, 0, momenta);

  if (true) {
    DiagrammaticMonteCarlo DMC(
      externalMomentum,
      maxLength,
      alpha,
      mu,
      numIterations,
      param,
      argv
    );
  } else {
    MonteCarlo{
      externalMomentum, 
      alpha,
      mu,
      numIterations,
      0,
      maxLength,
      argv
    };
  }


  printf("[Program finished @ %i in %.2fs]\n", worldRank, (double)(clock() - tStart)/CLOCKS_PER_SEC);


  MPI_Finalize();
  return 0;
}





