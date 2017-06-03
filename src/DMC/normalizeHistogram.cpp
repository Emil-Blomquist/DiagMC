#include "DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::normalizedHistogram (ArrayXXd& hist) {
  // give correct dimensions
  hist = ArrayXXd::Zero(this->hist.rows(), this->hist.cols());

  // calculate scale factor
  ArrayXd ps = ArrayXd::LinSpaced(this->hist.rows(), 0, this->maxMomenta - this->dp) + 0.5*this->dp;
  ArrayXd ts = ArrayXd::LinSpaced(this->hist.cols(), 0, this->maxLength - this->dt) + 0.5*this->dt;

  double sumG0 = 0;
  for (unsigned int i = 0; i != ps.size(); i++) {
    for (unsigned int j = 0; j != ts.size(); j++) {
      sumG0 += exp((this->mu - 0.5*ps[i]*ps[i])*ts[j]);

      // add dG contribution if we are going bold
      if (this->bold && this->boldIteration > 0) {
        sumG0 += this->dG(i, j);
      }
    }
  }

  double scaleFactor = sumG0/this->N0;

  // histogram corresponding to higher order diagrams
  hist = this->hist.cast<double>()*scaleFactor;
}

void DiagrammaticMonteCarlo::normalizedHistogram (ArrayXd& hist) {
  // give correct dimensions
  hist = ArrayXd::Zero(this->bins.size());

  // calcualte scale factor
  ArrayXd ts = ArrayXd::LinSpaced(this->bins.size(), 0, this->maxLength - this->dt) + 0.5*this->dt;

  double
    sumG0 = 0,
    p = this->initialExternalMomentum.norm();

  for (unsigned int j = 0; j != ts.size(); j++) {
    sumG0 += exp((this->mu - 0.5*p*p)*ts[j]);

    // add dG contribution if we are going bold
    if (this->bold && this->boldIteration > 0) {
      unsigned int pi = p/this->dp;
      sumG0 += this->dG(pi, j);
    }
  }

  double scaleFactor = sumG0/this->N0;

  // histogram corresponding to higher order diagrams
  hist = this->bins.cast<double>()*scaleFactor;
}