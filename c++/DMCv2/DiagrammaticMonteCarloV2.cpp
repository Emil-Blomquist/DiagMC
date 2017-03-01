#include "DiagrammaticMonteCarloV2.h"

DiagrammaticMonteCarloV2::DiagrammaticMonteCarloV2 (
  Vector3d P,
  double maxLength,
  double alpha,
  double mu,
  unsigned int numIterations,
  unsigned int maxOrder,
  unsigned int numBins,
  double param
) : FD{P, maxLength/2*2, alpha, mu} {
  // seed random generator
  unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
  this->mt.seed(seed1);

  // will print overflow exceptions etc.
  this->debug = true;
  // will print acceptance ratios etc. (debug must be true)
  this->loud = false;

  this->mu = mu;
  this->maxLength = maxLength;
  this->numIterations = numIterations;
  this->maxOrder = maxOrder;
  this->numBins = numBins;
  this->param = param;

  this->lastKeyMargin = 0.1;

  // store time at which the calculation began
  time_t rawtime;
  time (&rawtime);
  this->timeinfo = localtime(&rawtime);


  this->write2file();
  this->run();
}


void DiagrammaticMonteCarloV2::run () {

  // to reach start connfiguration
  const unsigned int untilStart = 10000;

  // to save data under the process
  const unsigned int saveAfter = 10*1000000;

  // Display display(&this->FD);

  // vector of pointers to member function of Phonon
  vector<double (DiagrammaticMonteCarloV2::*)(double)> updateMethods = {
    &DiagrammaticMonteCarloV2::shiftVertexPosition,
    &DiagrammaticMonteCarloV2::swapPhononConnections,
    &DiagrammaticMonteCarloV2::changeInternalPhononMomentumDirection,
    &DiagrammaticMonteCarloV2::changeInternalPhononMomentumMagnitude,
    &DiagrammaticMonteCarloV2::raiseOrder,
    &DiagrammaticMonteCarloV2::lowerOrder
    // &DiagrammaticMonteCarloV2::updateDiagramLength
  };





  // bins and keys for counting
  this->keys = vector<double>(this->numBins, 0);
  this->bins = vector<double>(this->numBins, 0);

  // fill keys
  double dt = this->maxLength/this->numBins;
  this->keys[this->numBins - 1] += this->lastKeyMargin; // so we never get out of bounds
  for (int i = 0; i != this->numBins; ++i) {
    this->keys[i] += (i + 1)*dt;
  }


  ///
  //
  ///

  unsigned int maxO = 0;



  for (int i = 0; i < this->numIterations + untilStart; ++i) {

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



    if (maxO < this->FD.Ds.size()) {
      maxO = this->FD.Ds.size();
      cout << maxO << endl;
    }

    if (i >= untilStart) {

      if ((i - untilStart)%saveAfter == saveAfter - 1) {
        // save temporary result
        this->write2file(i - untilStart + 1);
      }

      vector<double>::iterator itr = upper_bound(keys.begin(), keys.end(), this->FD.length);
      bins[itr - keys.begin()]++;
    }
  }

  // save final result
  this->write2file();
}


void DiagrammaticMonteCarloV2::write2file (const unsigned int iterationNum) {
  // create date and time string
  char buffer[80];
  strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", this->timeinfo);
  string dateAndTimeString(buffer);
  

  // create file name
  stringstream stream;
  stream << fixed << setprecision(3)
         << "p=" << this->FD.externalMomentum.norm()
         << " tmax=" << this->maxLength
         << " a=" << this->FD.couplingConstant
         << " mu=" << this->FD.chemicalPotential
         << " N=" << numIterations
         << " date=" << dateAndTimeString
         << " unique=" << param;
  string fileName = stream.str();

  // write to file
  ofstream myfile;
  myfile.open("../data/" + fileName + ".txt", ios_base::app);
  myfile << "-------- " + fileName;
  if (iterationNum) {
    myfile << " Ntemp=" << iterationNum;
  }
  myfile << " --------" << endl;
  myfile.close();

  //
  // TODO: we need to multiply with a factor (bins[0] new 1!)
  //



  // data to write to file
  myfile.open("../data/" + fileName + ".txt", ios_base::app);
  // for (auto key : this->keys) {
  //   if (key == this->keys[this->numBins - 1]) {
  //     myfile << fixed << setprecision(7) << key - this->lastKeyMargin <<  "\n";
  //   } else {
  //     myfile << fixed << setprecision(7) << key <<  " ";
  //   }
  // }
  // for (auto bin : this->bins) {
  //   if (bin == this->bins[this->numBins - 1]) {
  //     myfile << fixed << setprecision(7) << bin/this->bins[0] << "\n";
  //   } else {
  //     myfile << fixed << setprecision(7) << bin/this->bins[0] <<  " ";
  //   }
  // }

  if (! this->keys.empty()) {
    for (int i = 0; i != this->numBins; ++i) {
      if (i == this->numBins - 1) {
        myfile << fixed << setprecision(7) << this->keys[i] - this->lastKeyMargin <<  "\n";
      } else {
        myfile << fixed << setprecision(7) << this->keys[i] <<  " ";
      }
    }
    for (int i = 0; i != this->numBins; ++i) {
      if (i == this->numBins - 1) {
        myfile << fixed << setprecision(7) << this->bins[i]/this->bins[0] << "\n";
      } else {
        myfile << fixed << setprecision(7) << this->bins[i]/this->bins[0] <<  " ";
      }
    }
  }



  myfile.close();
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
  Vector3d tempVector, Ep, Eo1, Eo2, Qp, Qo, Q;
  
  Ep = P0.normalized();
  if (isfinite(Ep[0])) {
    // find a vector orthogonal to Ep
    tempVector = Ep.cross(Ep + Vector3d(1, 0, 0));
    Eo1 = tempVector.normalized();
    // cout << "TEMP " << Eo1.transpose() << endl;
    if ( ! isfinite(Eo1[0])) {
      Eo1  = Ep.cross(Ep + Vector3d(0, 1, 0)).normalized();
    }

    // span rest of R^3
    Eo2 = Ep.cross(Eo1);
  } else {
    // if P0 ≈ (0,0,0) we use that theta is the angle against the z-axis
    Ep = Vector3d(0, 0, 1);
    Eo1 = Vector3d(1, 0, 0);
    Eo2 = Vector3d(0, 1, 0);
  }

  // calculate the parallel and the orthogonal component of Q
  Qp = q*cos(theta)*Ep;
  Qo = (Eo1*cos(phi) + Eo2*sin(phi)) * q*sin(theta);
  Q = Qp + Qo;

  // to find overflow cause
  if (! isfinite(Q[0]) || ! isfinite(Q[1]) || ! isfinite(Q[2])) {
    cout << "-----------" << endl
         << "Overflow at DiagrammaticMonteCarloV2::calculateQ" << endl
         << "Q=" << Q.transpose() << endl
         << "P0= " << P0.transpose() << endl
         << "Ep= " << Ep.transpose() << endl
         << "Eo1= " << Eo1.transpose() << endl
         << "Eo2= " << Eo2.transpose() << endl
         << "q= " << q << endl
         << "theta= " << theta << endl
         << "phi= " << phi << endl
         << "-----------" << endl;
  }

  return Q;
}