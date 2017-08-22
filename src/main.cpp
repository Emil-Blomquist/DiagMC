#include <mpi.h> // MPI
#include <iostream>
#include <fstream>
#include <time.h> // execution time

#include "DMC/DiagrammaticMonteCarlo.h"
#include "MC/MonteCarlo.h"

using namespace std;

int main (int argc, char **argv) {
  clock_t tStart = clock();

  string pathToConfigFile = "";

  // input parameters
  int i = 1;
  while (i < argc) {
    if (strcmp(argv[i], "-config") == 0) {
      pathToConfigFile = argv[i + 1];
      i += 2;
    } else {
      i++;
    }
  }

  if (pathToConfigFile == "") {
    cout << ">>>Unable to run DMC since no config file given<<<" << endl;
    exit(EXIT_FAILURE);
  } else {
    string path = argv[0]; // program full path + name of binary file
    path.erase(path.find_last_of('/') + 1); // remove name of binary file

    ifstream configFile(path + "../" + pathToConfigFile);
    if ( ! configFile.good()) {
      cout << ">>>Unable to find the config file<<<" << endl;
      exit(EXIT_FAILURE);
    }
  }


  int worldRank, worldSize;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  if (true) {
    DiagrammaticMonteCarlo DMC(argv, pathToConfigFile);
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