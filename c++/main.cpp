#include <iostream>
#include "DMC/DiagrammaticMonteCarlo.h"

int main() {

  vector<double> times(250, 0);
  for(int i = 0; i != times.size(); i++) {
    times[i] = 0.01 + i*5.0/(times.size() - 1);
  }

  Vector3d externalMomentum(1, 0, 0);
  double 
    alpha = 2,
    mu = -2.2;

  int N = 10000000;


  for(int i = 79; i != times.size(); i++) {

    DiagrammaticMonteCarlo DMC(externalMomentum, times[i], alpha, mu);


    double g0 = DMC.FD();

    auto orderWeights = DMC.run(N, 7);

    // precision and shit
    cout << std::setprecision(7);
    cout << fixed;

    cout << times[i] << " ";
    for (auto order : orderWeights) {
      cout << g0*order/orderWeights[0] <<  " ";
    }
    cout << g0/orderWeights[0] << endl;
  }
}