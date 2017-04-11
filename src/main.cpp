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

  unsigned long int numIterations = 4000000000;

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

  // 50 bins per unit time
  const unsigned int numBins = 50*maxLength;

  Vector3d externalMomentum(0, 0, momenta);


  // Todo:
  // 1. Dyson equation
  //    1. Hash table of electron momenta
  //    2. determine proper diagrams
  // 2. function space
  // 3. skeleton diagrams


  // unordered_multimap<double, double> myMap;

  // myMap.insert({123.0, 1.0});
  // myMap.insert({123.0, 1.1});
  // myMap.insert({123.0, 1.2});
  // myMap.insert({123.1, 1.2});

  // // for (auto& x : myMap)
  // //   std::cout << x.first << ": " << x.second << std::endl;


  // auto range = myMap.equal_range(123.0);
  // for_each (
  //   range.first,
  //   range.second,
  //   [](pair<const double, double> x) {
  //     std::cout << x.second << endl;
  //   }
  // );


  DiagrammaticMonteCarlo DMC(
    externalMomentum,
    maxLength,
    alpha,
    mu,
    numIterations,
    numBins,
    param,
    argv
  );

  // MonteCarlo{
  //   externalMomentum, 
  //   alpha,
  //   mu,
  //   numIterations,
  //   maxLength,
  //   50,
  //   argv
  // };



  printf("[Finished in %.2fs]\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  return 0;
}


// #include <mpi.h> // MPI
// int myrank, nprocs;
// MPI_Init(NULL, NULL);
// MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
// MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

// ... code ...

// MPI_Finalize();