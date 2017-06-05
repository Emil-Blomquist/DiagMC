#include "DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::sumHistograms (
  Array<unsigned long long int, Dynamic, Dynamic>& totHist,
  unsigned long long int& totN0,
  unsigned long long int& totCurrItr
) {
  if (this->worldRank == 0) {
    // create a large enough buffer
    totHist = Array<unsigned long long int, Dynamic, Dynamic>::Zero(this->Np, this->Nt);
  }

  // collect total number of iterations
  MPI_Reduce(
    &this->currItr,
    &totCurrItr,
    1,
    MPI::UNSIGNED_LONG_LONG,
    MPI_SUM,
    0,
    MPI_COMM_WORLD);

  // collect counts for orders > 0
  MPI_Reduce(
    this->hist.data(),
    totHist.data(),
    this->hist.size(),
    MPI::UNSIGNED_LONG_LONG,
    MPI_SUM,
    0,
    MPI_COMM_WORLD);

  // collect counts for order 0
  MPI_Reduce(
    &this->N0,
    &totN0,
    1,
    MPI::UNSIGNED_LONG_LONG,
    MPI_SUM,
    0,
    MPI_COMM_WORLD);
}