#include "DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::normalizeHistogram ( 
  Array<unsigned long long int, Dynamic, Dynamic>& hist,
  unsigned long long int& N0,
  ArrayXXd& normalizedHist
) {
  // give correct dimensions
  normalizedHist = ArrayXXd::Zero(hist.rows(), hist.cols());

  // in the case that we have imported a G and we still have a fixed external momentum
  unsigned int Np = (this->fixedExternalMomentum ? 1 : this->Np);
  
  // calculate scale factor
  double sumG0 = 0;
  for (unsigned int i = 0; i != Np; i++) {
    double p = (this->fixedExternalMomentum ? initialExternalMomentum.norm() : (i + 0.5)*this->dp);
    for (unsigned int j = 0; j != this->Nt; j++) {
      double t = (j + 0.5)*this->dt;
      sumG0 += exp((this->mu - 0.5*p*p)*t);
    }
  }

  if (dG.size() > 0) {
    // add dG contribution
    if (this->fixedExternalMomentum) {
      unsigned int pi = min(initialExternalMomentum.norm()/this->dp, this->Np - 1.0);
      sumG0 += this->dG.row(pi).sum();
    } else {
      sumG0 += this->dG.sum();
    }
  }

  double scaleFactor = sumG0/N0;

  // histogram corresponding to higher order diagrams
  normalizedHist = hist.cast<double>()*scaleFactor;

  // add MC contribution from S1 diagram
  if (this->minDiagramOrder <= 1 && ! this->externalLegs) {
    for (unsigned int i = 0; i != this->S1mc.rows(); i++) {
      for (unsigned int j = 0; j != this->S1mc.cols(); j++) {
        normalizedHist(i, j) += this->S1mc(i, j);
      }
    }
  }
}