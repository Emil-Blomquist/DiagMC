#include "DiagrammaticMonteCarloV1.h"

DiagrammaticMonteCarloV1::DiagrammaticMonteCarloV1 (
  Vector3d P,
  double length,
  double alpha,
  double mu,
  double param
) : FD{P, length, alpha, mu} {
  // seed random generator
  unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
  this->mt.seed(seed1);

  // will print overflow exceptions etc.
  this->debug = true;
  // will print acceptance ratios etc. (debug must be true)
  this->loud = false;

  this->param = param;

  this->mu = mu;

  // set this outside
  this->maxLength = 5;
}


vector<double> DiagrammaticMonteCarloV1::run (int numIterations, int maxOrder) {
  this->maxOrder = maxOrder;

  // to reach start connfiguration
  int untilStart = 0;

  // Display display(&this->FD);

  // vector of pointers to member function of Phonon
  vector<double (DiagrammaticMonteCarloV1::*)(double)> updateMethods = {
    // &DiagrammaticMonteCarloV1::shiftVertexPosition,
    // &DiagrammaticMonteCarloV1::swapPhononConnections,
    // &DiagrammaticMonteCarloV1::changeInternalPhononMomentumDirection,
    // &DiagrammaticMonteCarloV1::changeInternalPhononMomentumMagnitude,
    // &DiagrammaticMonteCarloV1::raiseOrder,
    // &DiagrammaticMonteCarloV1::lowerOrder,
    &DiagrammaticMonteCarloV1::updateDiagramLength
  };

  // bins for counting
  vector<double> bins(this->maxOrder + 1, 0);

  for (int i = 0; i < numIterations + untilStart; ++i) {
    // choose update operation on random
    auto updateMethod = updateMethods[this->Uint(0, updateMethods.size() - 1)];

    // save configuration
    this->FD.save();

    // update diagram
    double a = (this->*updateMethod)(this->param);

    if (a < this->Udouble(0, 1)) {
      // rejected update
      this->FD.revert();
      if (updateMethod == &DiagrammaticMonteCarloV1::raiseOrder && this->FD.Ds.size() < this->maxOrder) {
        // unlink must happen here to all elements being removed
        this->phonon2beRemoved->unlink();
        this->vertices2beRemoved[0]->unlink();
        this->vertices2beRemoved[1]->unlink();
        this->electrons2beRemoved[0]->unlink();
        this->electrons2beRemoved[1]->unlink();
      }
    } else if (updateMethod == &DiagrammaticMonteCarloV1::lowerOrder) {
      // unlink must happen here to all elements being removed
      this->phonon2beRemoved->unlink();
      this->vertices2beRemoved[0]->unlink();
      this->vertices2beRemoved[1]->unlink();
      this->electrons2beRemoved[0]->unlink();
      this->electrons2beRemoved[1]->unlink();
    }

    // display.render();

    if (i >= untilStart) {
      bins[this->FD.Ds.size()]++;
    }
  }

  for(vector<double>::size_type i = 0; i != bins.size(); i++) {
    bins[i] /= numIterations;
  }

  return bins;
}














double DiagrammaticMonteCarloV1::Udouble (double from, double to) {
  uniform_real_distribution<double> distribution(from, to);
  return distribution(this->mt);
}

int DiagrammaticMonteCarloV1::Uint (int fromIncluded, int toIncluded) {
  uniform_int_distribution<int> distribution(fromIncluded, toIncluded);
  return distribution(this->mt);
}

Vector3d DiagrammaticMonteCarloV1::calculateP0 (shared_ptr<Phonon> d) {
  double dt = d->end->position - d->start->position;

  Vector3d P0(0, 0, 0);

  // loop through electrons between phonon propagator ends
  shared_ptr<Electron> g = d->start->G[1];
  while (g->start != d->end) {
    P0 += g->momentum*(g->end->position - g->start->position);

    g = g->end->G[1];
  }

  // add own momentum as well
  P0 = P0/dt + d->momentum;

  return P0;
}

Vector3d DiagrammaticMonteCarloV1::calculateQ (Vector3d P0, double q, double theta, double phi) {
  Vector3d tempVector, Ep, Eo1, Eo2, Qp, Qo;
  
  Ep = P0.normalized();
  if (::isnan(Ep[0])) { // ::isnan to prevent compiling error
    // if P0 ≈ (0,0,0) we use that theta is the angle against the z-axis
    Ep = Vector3d(0, 0, 1);
    Eo1 = Vector3d(1, 0, 0);
    Eo2 = Vector3d(0, 1, 0);
  } else {
    // find a vector orthogonal to Ep
    tempVector = Ep.cross(Ep + Vector3d(1, 0, 0));
    Eo1 = tempVector.normalized();
    // cout << "TEMP " << Eo1.transpose() << endl;
    if (::isnan(Eo1[0])) {
      Eo1  = Ep.cross(Ep + Vector3d(0, 1, 0)).normalized();
    }

    // span rest of R^3
    Eo2 = Ep.cross(Eo1);
  }

  // calculate the parallel and the orthogonal component of Q
  Qp = q*cos(theta)*Ep;
  Qo = (Eo1*cos(phi) + Eo2*sin(phi)) * q*sin(theta);

  return Qp + Qo;
}