#include "DiagrammaticMonteCarloV2.h"

DiagrammaticMonteCarloV2::DiagrammaticMonteCarloV2 (
  Vector3d P,
  double maxLength,
  double alpha,
  double mu,
  double param
) : FD{P, maxLength/2, alpha, mu} {
  // seed random generator
  unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
  this->mt.seed(seed1);

  // will print overflow exceptions etc.
  this->debug = true;
  // will print acceptance ratios etc. (debug must be true)
  this->loud = false;

  this->mu = mu;
  this->maxLength = maxLength;

  this->param = param;

  this->lastKeyMargin = 0.1;
}


tuple<vector<double>, vector<double>> DiagrammaticMonteCarloV2::run (
  unsigned int numIterations,
  unsigned int maxOrder,
  unsigned int numBins
) {
  this->maxOrder = maxOrder;

  // to reach start connfiguration
  int untilStart = 0;

  // Display display(&this->FD);

  // vector of pointers to member function of Phonon
  vector<double (DiagrammaticMonteCarloV2::*)(double)> updateMethods = {
    &DiagrammaticMonteCarloV2::shiftVertexPosition,
    &DiagrammaticMonteCarloV2::swapPhononConnections,
    &DiagrammaticMonteCarloV2::changeInternalPhononMomentumDirection,
    &DiagrammaticMonteCarloV2::changeInternalPhononMomentumMagnitude,
    &DiagrammaticMonteCarloV2::raiseOrder,
    &DiagrammaticMonteCarloV2::lowerOrder,
    &DiagrammaticMonteCarloV2::updateDiagramLength
  };





  // bins and keys for counting
  vector<double>
    keys(numBins, 0),
    bins(numBins, 0);

  // fill keys
  double dt = this->maxLength/numBins;
  keys[numBins - 1] += this->lastKeyMargin; // so we never get out of bounds
  for (int i = 0; i != keys.size(); ++i) {
    keys[i] += (i + 1)*dt;
  }






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
      if (updateMethod == &DiagrammaticMonteCarloV2::raiseOrder && this->FD.Ds.size() < this->maxOrder) {
        // unlink must happen here to all elements being removed
        this->phonon2beRemoved->unlink();
        this->vertices2beRemoved[0]->unlink();
        this->vertices2beRemoved[1]->unlink();
        this->electrons2beRemoved[0]->unlink();
        this->electrons2beRemoved[1]->unlink();
      }
    } else if (updateMethod == &DiagrammaticMonteCarloV2::lowerOrder) {
      // unlink must happen here to all elements being removed
      this->phonon2beRemoved->unlink();
      this->vertices2beRemoved[0]->unlink();
      this->vertices2beRemoved[1]->unlink();
      this->electrons2beRemoved[0]->unlink();
      this->electrons2beRemoved[1]->unlink();
    }

    // display.render();

    if (i >= untilStart) {
      // bins[(int)(this->FD.length/dt)]++;

      vector<double>::iterator itr = upper_bound(keys.begin(), keys.end(), this->FD.length);
      bins[itr - keys.begin()]++;
    }
  }

  //
  // TODO: we need to multiply with a factor (bins[0] new 1!)
  //

  // prepare result
  double firstBin = bins[0];
  for(vector<double>::size_type i = 0; i != bins.size(); i++) {
    bins[i] /= firstBin;
  }
  keys[numBins - 1] -= this->lastKeyMargin;

  return make_tuple(keys, bins);
}














double DiagrammaticMonteCarloV2::Udouble (double from, double to) {
  uniform_real_distribution<double> distribution(from, to);
  return distribution(this->mt);
}

int DiagrammaticMonteCarloV2::Uint (int fromIncluded, int toIncluded) {
  uniform_int_distribution<int> distribution(fromIncluded, toIncluded);
  return distribution(this->mt);
}

Vector3d DiagrammaticMonteCarloV2::calculateP0 (shared_ptr<Phonon> d) {
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

Vector3d DiagrammaticMonteCarloV2::calculateQ (Vector3d P0, double q, double theta, double phi) {
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