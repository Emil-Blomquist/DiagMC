#include <iostream>
#include <time.h> // execution time

#include <unordered_set>


#include "DMCv2/DiagrammaticMonteCarloV2.h"

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

  long int numIterations = 4000000000;

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

  // DiagrammaticMonteCarloV2 DMC(
  //   externalMomentum,
  //   maxLength,
  //   alpha,
  //   mu,
  //   numIterations,
  //   numBins,
  //   param,
  //   argv
  // );


  unordered_set<shared_ptr<int> > myset;

  shared_ptr<int>
    p1(new int{123}),
    p2(new int{321});

  myset.insert(p1);
  myset.insert(p2);

  cout << *p1 << " " << *p2 << endl;

  *p1 = 666;

  cout << *p1 << " " << *p2 << endl;


  for (const shared_ptr<int>& x: myset) std::cout << " " << *x;





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