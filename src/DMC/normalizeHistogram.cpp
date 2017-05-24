#include "DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::normalizedHistogram (Array<double, Dynamic, Dynamic>& hist) {
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

  // to remove the error due the combination of a singular diagram and a discretized time
  vector<double> singularityFix((int) round(this->maxLength/this->dt), 0);
  // if ( ! this->externalLegs && this->minDiagramOrder <= 1) {
  //   for (unsigned int i = 0; i != singularityFix.size(); ++i) {
  //     singularityFix[i] = this->alpha*(
  //       exp(-0.5*(2*i + 1)*this->dt)/sqrt(M_PI*0.5*(2*i + 1)*this->dt)
  //       - (erf(sqrt((i + 1)*this->dt)) - erf(sqrt(i*this->dt)))/this->dt
  //     );
  //   }
  // }

  // histogram corresponding to higher order diagrams
  for (unsigned int i = 0; i != this->hist.rows(); i++) {
    for (unsigned int j = 0; j != this->hist.cols(); j++) {
      hist(i, j) = abs(this->hist(i, j)*scaleFactor + singularityFix[j]);
    }
  }
}