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
      maxLength,
      50,
      argv
    };
  }

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


// 1D och 2D överenstämmer ej. 2D verkar inte heller överenstämma mot andra 2D...
//   1. kör 1D i lång tid på clustret och se om long time 1D överenstämmer med 2D 
//   2. kan ha blivit knas när ny update funktion implementerades
//   3. testa att köra 2D fast med mindre pMax
//   4. variera dp
// Liknande problem uppstod när vi introducerade irreducible diagrams
//   1. testa att bara beräkna den vanliga greens functionen?
//   2. irreducible + external legs, funkar det? -> om inte, irreducible funkar inte!
//   3. Om G0 verkar funka, skit it och gå vidare istället
// Kolla om updateTime ger a - 1 ≈10^-6 när man försöker sampla G0