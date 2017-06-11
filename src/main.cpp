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
    cout << endl;
    cout << "check the statistics of the first iteration of BOLD where dE = 0. Is the statistic as bad?" << endl;
    cout << endl;
  }







  // cout << endl
  //      << "-7) Calculate with MC using dE the first order self energy" << endl
  //      << "-6) Since we are only doing first order, use irreducible instead of skeleton diagram" << endl
  //      << "-5) Compare the self energy with the known one" << endl
  //      << "-4) Turn off all unimportant update function for when only using 1st order" << endl
  //      << "-3) If this is due to the first order diagram, that would be weird since we have succesfully done Dyson in 2D" << endl
  //      << "-2) Compare the G's of the first iteration running on the cluster" << endl
  //      << "-1) Instead of saving dG, save S" << endl
  //      << "0) Make sure that the files updates update their content" << endl
  //      << "0.5) Perhaps we have extremely poor distribution of parameters for bold algorithm? Maybe implement completely new set of update functions?" << endl
  //      << "1) Import G from file and calculate corresponding dE." << endl
  //      << "2) Check that using this dE, G is always obtained, no matter what order we allow." << endl
  //      << "3) Starting from dE = 0 and allowing only first order, check so that the G obtained from bold is larger then the one obtained from Dyson." << endl
  //      << "4) Perhaps has something to do with dG/G0 when we create dE?" << endl
  //      << endl;

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





