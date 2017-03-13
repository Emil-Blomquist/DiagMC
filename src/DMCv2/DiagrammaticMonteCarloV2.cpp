#include "DiagrammaticMonteCarloV2.h"

DiagrammaticMonteCarloV2::DiagrammaticMonteCarloV2 (
  Vector3d P,
  double maxLength,
  double alpha,
  double mu,
  unsigned long int numIterations,
  unsigned int numBins,
  double param,
  char **argv
) : FD(P, maxLength/2, alpha, mu) {
  // seed random generator
  unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
  this->mt.seed(seed1);

  // will print overflow exceptions etc. (code is half as fast with this activated)
  this->debug = false;
  // will print acceptance ratios etc. (debug must be true)
  this->loud = false;

  this->mu = mu;
  this->alpha = alpha;
  this->maxLength = maxLength;
  this->numIterations = numIterations;
  this->numBins = numBins;
  this->param = param;
  this->argv = argv;

  this->lastKeyMargin = 0.1;

  // number of times we visit the zeroth order diagram of the zeroth bin
  this->n00 = 0;

  // store time at which the calculation began
  time_t rawtime;
  time (&rawtime);
  this->timeinfo = localtime(&rawtime);


  this->write2file();
  this->run();
}


void DiagrammaticMonteCarloV2::run () {
  // to reach start connfiguration
  const unsigned int untilStart = 100000;

  // to save data under the process
  const unsigned int saveAfter = 50*1000000;

  // vector of pointers to member function of Phonon
  vector<void (DiagrammaticMonteCarloV2::*)(double)> updateMethods = {
    &DiagrammaticMonteCarloV2::shiftVertexPosition, // <- 2
    &DiagrammaticMonteCarloV2::swapPhononConnections, // <- 1
    &DiagrammaticMonteCarloV2::changeInternalPhononMomentumDirection, // <- 3
    &DiagrammaticMonteCarloV2::changeInternalPhononMomentumMagnitude, // <- 4
    &DiagrammaticMonteCarloV2::raiseOrder,
    &DiagrammaticMonteCarloV2::lowerOrder,
    &DiagrammaticMonteCarloV2::changeDiagramLength
  };

  // bins and keys for counting
  this->keys = vector<double>(this->numBins, 0);
  this->bins = vector<int>(this->numBins, 0);

  // fill keys
  double dt = this->maxLength/this->numBins;
  this->keys[this->numBins - 1] += this->lastKeyMargin; // so we never get out of bounds
  for (int i = 0; i != this->numBins; ++i) {
    this->keys[i] += (i + 1)*dt;
  }

  // to start at a random position
  for (unsigned int i = 0; i < untilStart; ++i) {
    // choose update operation on random
     auto updateMethod = updateMethods[this->Uint(0, updateMethods.size() - 1)];
    // update diagram
    (this->*updateMethod)(this->param);
  }


  // Display disp(&this->FD);
  // disp.render();


  // main loop
  for (unsigned long int i = 0; i < this->numIterations; ++i) {
    auto updateMethod = updateMethods[this->Uint(0, updateMethods.size() - 1)];
    (this->*updateMethod)(this->param);

    if (i%saveAfter == saveAfter - 1) {
      // save temporary result
      this->write2file(i + 1);
    }

    // bin diagram length
    vector<double>::iterator itr = upper_bound(keys.begin(), keys.end(), this->FD.length);
    int index = itr - keys.begin();
    bins[index]++;

    // if at first bin and at zeroth order, bin again
    if (index == 0 && this->FD.Ds.size() == 0) {
      this->n00++;
    }

  }

  // save final result
  this->write2file();
}


void DiagrammaticMonteCarloV2::write2file (const unsigned long int iterationNum) {
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

  // obtain path relative binary file
  string path = this->argv[0]; // program full path + name of binary file
  path.erase(path.find_last_of('/') + 1); // remove name of binary file

  // create data folder if not already existing
  system(("mkdir -p " + path + "../data").c_str());

  // write to file
  ofstream myfile;
  myfile.open(path + "../data/" + fileName + ".txt");
  myfile << "-------- " + fileName;
  if (iterationNum) {
    myfile << " Ntemp=" << iterationNum;
  }
  myfile << " --------" << endl;
  myfile.close();

  if (! this->keys.empty()) {
    // data to write to file
    myfile.open(path + "../data/" + fileName + ".txt", ios_base::app);

    for (int i = 0; i != this->numBins; ++i) {
      if (i == 0) {
        myfile << fixed << setprecision(7) << 0.5*this->keys[i] << " ";
      } else if (i == this->numBins - 1) {
        myfile << fixed << setprecision(7) << 0.5*(this->keys[i - 1] + this->keys[i] - this->lastKeyMargin) <<  "\n";
      } else {
        myfile << fixed << setprecision(7) << 0.5*(this->keys[i - 1] + this->keys[i]) <<  " ";
      }
    }

    double
      g0Ofdt = exp(-0.5*this->keys[0]*(0.5*this->FD.externalMomentum.squaredNorm() - this->mu)),
      gOfdt = (double) this->bins[0]/this->n00 * g0Ofdt;

    for (int i = 0; i != this->numBins; ++i) {
      if (i == this->numBins - 1) {
        myfile << fixed << setprecision(7) << (double) this->bins[i]/this->bins[0] * gOfdt;
      } else {
        myfile << fixed << setprecision(7) << (double) this->bins[i]/this->bins[0] * gOfdt <<  " ";
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

Vector3d DiagrammaticMonteCarloV2::calculateMeanP (shared_ptr<Vertex> v1, shared_ptr<Vertex> v2) {

  Vector3d meanP(0, 0, 0);

  if (v1 != v2) {
    if (v1->G[1]->end == v2) {
      // if dt = 0
      meanP += v1->G[1]->momentum;

      if ( ! isfinite(meanP[0]) || ! isfinite(meanP[1]) || ! isfinite(meanP[2])) {
        cout << "-----------" << endl
             << "Overflow at DiagrammaticMonteCarloV2::calculateMeanP" << endl
             << "meanP= " << meanP.transpose() << endl
             << "t1= " << v1->position << endl
             << "t2= " << v2->position << endl
             << "-----------" << endl;
      }
    } else {
      double dt = v2->position - v1->position;

      // loop through electrons between phonon propagator ends
      shared_ptr<Electron> g = v1->G[1];
      while (g->start != v2) {
        meanP += g->momentum*(g->end->position - g->start->position);

        g = g->end->G[1];
      }

      if (dt == 0 || ! isfinite(meanP[0]) || ! isfinite(meanP[1]) || ! isfinite(meanP[2])) {
        cout << "-----------" << endl
             << "Overflow at DiagrammaticMonteCarloV2::calculateMeanP" << endl
             << "dt=" << dt << endl
             << "meanP= " << meanP.transpose() << endl
             << "t1= " << v1->position << endl
             << "t2= " << v2->position << endl
             << "-----------" << endl;
      }

      meanP /= dt;
    }
  }

  return meanP;
}

Vector3d DiagrammaticMonteCarloV2::calculateP0 (shared_ptr<Phonon> d) {
  Vector3d meanP = this->calculateMeanP(d->start, d->end);

  return meanP + d->momentum;
}

Vector3d DiagrammaticMonteCarloV2::calculateQ (Vector3d P0, double q, double theta, double phi) {
  Vector3d tempVector, Ep, Eo1, Eo2, Qp, Qo, Q;
  
  Ep = P0.normalized();
  if (isfinite(Ep[0])) {
    // find a vector orthogonal to Ep
    tempVector = Ep.cross(Ep + Vector3d(1, 0, 0));
    Eo1 = tempVector.normalized();

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