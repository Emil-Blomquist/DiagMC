#include "DiagrammaticMonteCarlo.h"

double DiagrammaticMonteCarlo::dEOf (shared_ptr<Electron> g) {
  unsigned int
    pi = min(g->p/this->dp, this->dE.rows() - 1.0),
    ti = (g->end->position - g->start->position)/this->dt;

  return this->dE(pi, ti);
}

double DiagrammaticMonteCarlo::dEOf (double p, double t) {
  unsigned int
    pi = min(p/this->dp, this->dE.rows() - 1.0),
    ti = t/this->dt;

  return this->dE(pi, ti);
}

double DiagrammaticMonteCarlo::additionalPhase (shared_ptr<Electron> g) {
  double t = g->end->position - g->start->position;

  unsigned int
    pi = min(g->p/this->dp, this->dE.rows() - 1.0),
    ti = t/this->dt;

  return -this->dE(pi, ti)*t;
}

double DiagrammaticMonteCarlo::additionalPhase (Vector3d P, double t) {
  double p = P.norm();

  unsigned int
    pi = min(p/this->dp, this->dE.rows() - 1.0),
    ti = t/this->dt;

  return -this->dE(pi, ti)*t;
}

double DiagrammaticMonteCarlo::additionalPhase (double p, double t) {
  unsigned int
    pi = min(p/this->dp, this->dE.rows() - 1.0),
    ti = t/this->dt;

  return -this->dE(pi, ti)*t;
}

void DiagrammaticMonteCarlo::calculateEnergyDiff () {
  for (unsigned int i = 0; i != this->dE.rows(); i++) {
    double
      p = (i + 0.5)*this->dp,
      E0 = 0.5*p*p;

    for (unsigned int j = 0; j != this->dE.cols(); j++) {
      double
        t = (j + 0.5)*this->dt,
        G0 = exp((this->mu - E0)*t),
        logOfG = (this->mu - E0)*t + (isfinite(this->dG(i, j)/G0) ? log(1 + abs(this->dG(i, j)/G0)) : 0);

      this->dE(i, j) = this->mu - E0 - logOfG/t;

      if ( ! isfinite(this->dE(i, j))) {
        cout << p << " " << t << fixed << setprecision(16) << ": " << log(G0 + this->dG) << " vs "
             << (this->mu - 0.5*p*p)*t + log(1 + this->dG/G0) << " but " << this->dG/G0 << endl;
      }
    }
  }

  // // interpolate dE
  // for (unsigned int ii = 0; ii != this->dE.rows(); ii++) {
  //   double p = (ii + 0.5)*this->dp * this->dG.rows()/this->dE.rows();
  //   for (unsigned int jj = 0; jj != this->dE.cols(); jj++) {
  //     double t = (jj + 0.5)*this->dt * this->dG.cols()/this->dE.cols();

  //     double
  //       i = p/this->dp,
  //       j = t/this->dt;

  //     unsigned int
  //       iu = round(i),
  //       ju = round(j);

  //     if (iu < 1) {
  //       if (ju < 1) {
  //         // 00 corner
  //         this->dE(ii, jj) = dE_low(0, 0);
  //       } else if (ju > dE_low.cols() - 1) {
  //         // 01 corner
  //         this->dE(ii, jj) = dE_low(0, dE_low.cols() - 1);
  //       } else {
  //         // smooth using the first row only
  //         this->dE(ii, jj) = dE_low(0, ju - 1)*(0.5 - j + ju) + dE_low(0, ju)*(0.5 + j - ju);
  //       }
  //     } else if (iu > dE_low.rows() - 1) {
  //       if (ju < 1) {
  //         // 10 corner
  //         this->dE(ii, jj) = dE_low(dE_low.rows() - 1, 0);
  //       } else if (ju > dE_low.cols() - 1) {
  //         // 11 corner
  //         this->dE(ii, jj) = dE_low(dE_low.rows() - 1, dE_low.cols() - 1);
  //       } else {
  //         // smooth using the last row only
  //         this->dE(ii, jj) = dE_low(dE_low.rows() - 1, ju - 1)*(0.5 - j + ju) + dE_low(dE_low.rows() - 1, ju)*(0.5 + j - ju);
  //       }
  //     } else {
  //       if (ju < 1) {
  //         // smooth using the first column only
  //         this->dE(ii, jj) = dE_low(iu - 1, 0)*(0.5 - i + iu) + dE_low(iu, 0)*(0.5 + i - iu);
  //       } else if (ju > dE_low.cols() - 1) {
  //         // smooth using the last column only
  //         this->dE(ii, jj) = dE_low(iu - 1, dE_low.cols() - 1)*(0.5 - i + iu) + dE_low(iu, dE_low.cols() - 1)*(0.5 + i - iu);
  //       } else {
  //         // smooth using the nearest four values
  //         this->dE(ii, jj) = dE_low(iu - 1, ju - 1)*(0.5 - i + iu)*(0.5 - j + ju)
  //                        + dE_low(iu - 1, ju)*(0.5 - i + iu)*(0.5 + j - ju)
  //                        + dE_low(iu, ju - 1)*(0.5 + i - iu)*(0.5 - j + ju)
  //                        + dE_low(iu, ju)*(0.5 + i - iu)*(0.5 + j - ju);
  //       }
  //     }
  //   }
  // }
}