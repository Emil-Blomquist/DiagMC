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
    p = 0,
    alpha = 1,
    mu = -1.1,
    maxLength = 40,
    maxMomenta = 1,
    dt = 0.02,
    dp = 0.02,
    numSecsDMC = 10*60,
    numSecsMC = 1*60;
  
  unsigned int numTempSaves = 9;


  // input parameters
  int i = 1;
  while (i < argc) {
    if (strcmp(argv[i], "-p") == 0) {
      p = stod(argv[i + 1]);
      i += 2;
    } else if (strcmp(argv[i], "-a") == 0) {
      alpha = stod(argv[i + 1]);
      i += 2;
    } else if (strcmp(argv[i], "-mu") == 0) {
      mu = stod(argv[i + 1]);
      i += 2;
    } else if (strcmp(argv[i], "-max_t") == 0) {
      maxLength = stod(argv[i + 1]);
      i += 2;
    } else if (strcmp(argv[i], "-max_p") == 0) {
      maxMomenta = stod(argv[i + 1]);
      i += 2;
    }  else if (strcmp(argv[i], "-dt") == 0) {
      dt = stod(argv[i + 1]);
      i += 2;
    }  else if (strcmp(argv[i], "-dp") == 0) {
      dp = stod(argv[i + 1]);
      i += 2;
    }  else if (strcmp(argv[i], "-num_secs_dmc") == 0) {
      numSecsDMC = stod(argv[i + 1]);
      i += 2;
    }  else if (strcmp(argv[i], "-num_secs_mc") == 0) {
      numSecsMC = stod(argv[i + 1]);
      i += 2;
    }  else if (strcmp(argv[i], "-num_temp_saves") == 0) {
      numTempSaves = stoul(argv[i + 1]);
      i += 2;
    } else {
      i++;
    }
  }

  int worldRank, worldSize;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  if (true) {
    DiagrammaticMonteCarlo DMC(
      p,
      alpha,
      mu,
      maxLength,
      maxMomenta,
      dt,
      dp,
      numSecsDMC,
      numSecsMC,
      numTempSaves,
      argv
    );
  }
  // else {
  //   MonteCarlo{
  //     externalMomentum, 
  //     alpha,
  //     mu,
  //     numIterations,
  //     0,
  //     maxLength,
  //     argv
  //   };
  // }


  printf("[Program finished @ %i in %.2fs]\n", worldRank, (double)(clock() - tStart)/CLOCKS_PER_SEC);


  MPI_Finalize();
  return 0;
}