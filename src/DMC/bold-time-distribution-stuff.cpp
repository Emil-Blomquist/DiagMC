#include "DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::calculateLambdas (ArrayXXd& dE, ArrayXd& lambdas) {
  // create large enough buffer
  lambdas = ArrayXd::Zero(dE.rows());

  // times to be used in the fitting process
  ArrayXd t = ArrayXd::LinSpaced(this->Nt, 0, (this->Nt - 1)*this->dt) + 0.5*this->dt;

  // from each row, extract a decay
  for (unsigned int i = 0; i != dE.rows(); i++) {
    double p = (i + 0.5)*this->dp;

    // save row temporarily
    ArrayXd dEofP = dE.row(i);

    lambdas[i] = this->expFit(t, dEofP, p, this->mu);
  }
}

double DiagrammaticMonteCarlo::expFit(ArrayXd& t, ArrayXd& dE, double p, double mu) {
  // create G from t and dE
  ArrayXd G = exp((mu - 0.5*p*p - dE)*t);

  // threshold for elements to consider during the fitting
  const double threshold = 1.0/1000;

  // figure out which elements to keep
  unsigned int n = dE.size();
  for (unsigned int i = 1; i != G.size(); i++) {
    if (G[i] < threshold) {
      n = i + 1;
      break;
    }
  }

  // if we have less then two elements, we are unable to do a fit
  if (n == 1) {
    return 0.5*p*p - mu;
  }

  // keep the "n" first elements and create linear equation
  VectorXd y = log(G.head(n));

  MatrixXd A(n, 2);
  A.col(0) = MatrixXd::Constant(n, 1, 1);
  A.col(1) = -t.head(n);

  ArrayXd fittedParameters = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(y);

  double
    // C = exp(fittedParameters[0]),
    lambda = fittedParameters[1];

  return lambda;
}

double DiagrammaticMonteCarlo::lambdaOf (double p) {
  unsigned int pi = p/this->dp;

  // if p larger than maximum discretized one, return bare energy
  if (pi > this->dE.rows() - 1) {
    return 0.5*p*p - this->mu;
  } else {
    return this->lambdas[pi]; 
  }
}

double DiagrammaticMonteCarlo::lambdaOf (shared_ptr<Electron> g) {
  unsigned int pi = g->p/this->dp;

  // if p larger than maximum discretized one, return bare energy
  if (pi > this->dE.rows() - 1) {
    return 0.5*g->p*g->p - this->mu;
  } else {
    return this->lambdas[pi]; 
  }
}