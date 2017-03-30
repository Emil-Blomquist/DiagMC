#include <iostream>
#include <time.h> // execution time

#include "DMC/DiagrammaticMonteCarlo.h"

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
  unsigned int i = 1;
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

  //
  // 1. alla updates har olika sannolikhet och jämför olika resultat - gjort
  // 2. samma som 1. fast gör med flit så att en uppdatering blir felaktig och se om det skiljer - gjort
  // 3. bin size - gjort
  // 4. kolla olika riktingar på external momentum - gjort
  // 5. kolla momentum konservation efter varje uppdatering
  // 6. kolla de där konstiga potentialerna som krävdes när man skulle ändra tid enligt artikeln
  // 7. testa med hardcore mode fast med moment i z-rikting
  // 8. be om resultat
  // 9. jämför exakta resultat - gjort
  // 10. implementera bra raise/lower order - gjort
  //

  DiagrammaticMonteCarloV2 DMC(
    externalMomentum,
    maxLength,
    alpha,
    mu,
    numIterations,
    numBins,
    param,
    argv
  );

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