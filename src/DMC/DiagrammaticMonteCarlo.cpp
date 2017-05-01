#include "DiagrammaticMonteCarlo.h"

DiagrammaticMonteCarlo::DiagrammaticMonteCarlo (
  Vector3d P,
  double maxLength,
  double alpha,
  double mu,
  unsigned long int numIterations,
  // unsigned int numBins,
  double param,
  char **argv
) : FD(P, maxLength/2, alpha, mu) {
  // seed random generator
  unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
  this->mt.seed(seed1);

  // will print overflow exceptions etc. (code is half as fast with this activated)
  this->debug = false;
  // will print acceptance ratios etc. (debug must be true)
  this->loud = false;

  // if the diagram should have external legs or not
  this->externalLegs = false;
  // if we should only take into account irreducible diagrams
  this->reducibleDiagrams = false;
  // if we should let hte external momentum vary or not
  this->fixedExternalMomentum = false;
  // for when to bin the diagram
  this->minDiagramOrder = 1;
  // raise diagrm order will look at this (zero -> turned off)
  this->maxDiagramOrder = 0;

  this->mu = mu;
  this->alpha = alpha;
  this->numIterations = numIterations;
  this->param = this->Udouble(0, 1);
  this->argv = argv;

  // store time at which the calculation began
  time_t rawtime;
  time (&rawtime);
  this->timeinfo = localtime(&rawtime);

  // in order to create the histogram
  this->maxLength = maxLength;
  this->dt = 0.02;

  if ( ! this->fixedExternalMomentum) {
    this->maxMomenta = 0.001;
    this->dp = 0.001;

    const unsigned int
      Np = this->maxMomenta/this->dp,
      Nt = this->maxLength/this->dt;

    this->N0 = 0;
    this->hist = Array<unsigned long int, Dynamic, Dynamic>::Zero(Np, Nt);
  } else {
    this->bins = vector<unsigned long int>(this->maxLength/this->dt, 0);
    this->bins0 = vector<unsigned long int>(this->maxLength/this->dt, 0);
  }

  this->write2file();
  this->run();
}

void DiagrammaticMonteCarlo::run () {
  // to reach start connfiguration
  const unsigned int untilStart = 10000000;

  // to save data under the process
  const unsigned int saveAfter = 500*1000000;

  // specify the relative probability of choosing a specific update function
  multimap<unsigned int, void (DiagrammaticMonteCarlo::*)(double)> updateMethods = {
    {1, &DiagrammaticMonteCarlo::shiftVertexPosition},
    {1, &DiagrammaticMonteCarlo::swapPhononConnections},
    {1, &DiagrammaticMonteCarlo::changeInternalPhononMomentumDirection},
    {1, &DiagrammaticMonteCarlo::changeInternalPhononMomentumMagnitude},
    {1, &DiagrammaticMonteCarlo::raiseOrder}, // <- These two must have the same probability
    {1, &DiagrammaticMonteCarlo::lowerOrder}, // <-
    {1, &DiagrammaticMonteCarlo::changeDiagramLength},
    {1, &DiagrammaticMonteCarlo::changeDiagramLengthComplex},
    {(this->fixedExternalMomentum ? 0 : 1), &DiagrammaticMonteCarlo::changeExternalMomentumMagnitude}
  };

  // vector which is going to contain the specified quantity of update functions
  vector<void (DiagrammaticMonteCarlo::*)(double)> chooseUpdateMethod;

  // populate vector
  for (auto updateMethod = updateMethods.begin(); updateMethod != updateMethods.end(); updateMethod++) {
    for (unsigned int i = 0; i != updateMethod->first; i++) {
      chooseUpdateMethod.push_back(updateMethod->second);
    }
  }

  // to start at a random position
  for (unsigned int i = 0; i < untilStart; ++i) {
    // choose update operation on random
    auto updateMethod = chooseUpdateMethod[this->Uint(0, chooseUpdateMethod.size() - 1)];
    // update diagram
    (this->*updateMethod)(this->param);
  }

  // main loop
  for (unsigned long int i = 0; i < this->numIterations; ++i) {
    auto updateMethod = chooseUpdateMethod[this->Uint(0, chooseUpdateMethod.size() - 1)];
    (this->*updateMethod)(this->param);

    if (i%saveAfter == saveAfter - 1) {
      // save temporary result
      this->write2file(i + 1);
    }

    // bin diagrams of order higher than 0
    if (this->FD.Ds.size() >= this->minDiagramOrder) {
      if (this->reducibleDiagrams || this->FD.diagramIsIrreducible()) {
        if (this->fixedExternalMomentum) {
          unsigned int index = this->FD.length/this->dt;
          bins[index]++;
        } else {
          unsigned int
            pi = this->FD.externalMomentum/this->dp,
            ti = this->FD.length/this->dt;
          this->hist(pi, ti)++;
        }
      }
    }

    // if zeroth order, bin again
    if (this->FD.Ds.size() == 0) {
      if (this->fixedExternalMomentum) {
        unsigned int index = this->FD.length/this->dt;
        bins0[index]++;
      } else {
        this->N0++;
      }
    }
  }

  // save final result
  this->write2file();
}



void DiagrammaticMonteCarlo::write2file (const unsigned long int iterationNum) {
  // create date and time string
  char buffer[80];
  strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", this->timeinfo);
  string dateAndTimeString(buffer);

  // create file name
  stringstream stream;
  stream << fixed << setprecision(6) // match the precision of "to_string()"
         << " a=" << this->FD.couplingConstant
         << ( this->fixedExternalMomentum ? " p=" + to_string(this->FD.externalMomentum) : "" )
         << " mu=" << this->FD.chemicalPotential
         << " tmax=" << this->maxLength
         << " dt=" << this->dt
         << ( ! this->fixedExternalMomentum ? " pmax=" + to_string(this->maxMomenta) : "" )
         << ( ! this->fixedExternalMomentum ? " dp=" + to_string(this->dp) : "" )
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
  myfile.open(path + "../data/" + fileName.substr(1) + ".txt");
  myfile << "--------" + fileName;
  if (iterationNum) {
    myfile << " Ntemp=" << iterationNum;
  }
  myfile << " --------" << endl;
  myfile.close();


  // to remove the error due the combination of a singular diagram and a discretized time
  vector<double> singularityFix((int) this->maxLength/this->dt, 0);
  if ( ! this->externalLegs && this->minDiagramOrder <= 1) {
    for (unsigned int i = 0; i != this->maxLength/this->dt; ++i) {
      singularityFix[i] = this->alpha*(
        exp(-0.5*(2*i + 1)*this->dt)/sqrt(M_PI*0.5*(2*i + 1)*this->dt)
        - (erf(sqrt((i + 1)*this->dt)) - erf(sqrt(i*this->dt)))/this->dt
      );
    }
  }

  if ( ! this->fixedExternalMomentum && this->N0) {
    // open file
    myfile.open(path + "../data/" + fileName.substr(1) + ".txt", ios_base::app);

    // calculate scale factor
    ArrayXd ps = ArrayXd::LinSpaced(this->hist.rows(), 0, this->maxMomenta - this->dp) + 0.5*this->dp;
    ArrayXd ts = ArrayXd::LinSpaced(this->hist.cols(), 0, this->maxLength - this->dt) + 0.5*this->dt;

    double sumG0 = 0;
    for (unsigned int i = 0; i != ps.size(); i++) {
      for (unsigned int j = 0; j != ts.size(); j++) {
        sumG0 += exp((this->mu - 0.5*ps[i]*ps[i])*ts[j]);
      }
    }

    double scaleFactor = sumG0/this->N0;

    // histogram corresponding to higher order diagrams
    for (unsigned int i = 0; i != this->hist.rows(); i++) {
      for (unsigned int j = 0; j != this->hist.cols(); j++) {
        myfile << abs(this->hist(i, j)*scaleFactor + singularityFix[j]);
        if (j + 1 < this->hist.cols()) {
          myfile << " ";
        }
      }
      if (i + 1 < this->hist.rows()) {
        myfile << endl;
      }
    }
  } else if (this->fixedExternalMomentum && this->bins0[0] > 0) {
    // open file
    myfile.open(path + "../data/" + fileName.substr(1) + ".txt", ios_base::app);

    // times
    for (unsigned int i = 0; i != this->maxLength/this->dt; ++i) {
      myfile << fixed << setprecision(7) << (i + 0.5)*this->dt;
      if (i == this->maxLength/this->dt - 1) {
        myfile << "\n";
      } else {
        myfile << " ";
      }
    }

    // calculate scale factor
    unsigned int until = 0;
    while ((double) this->bins0[until]/this->bins0[0] > 0.01) {
      until++;
    }
    double Z = 0;
    for (unsigned int i = 0; i != until; ++i) {
      Z += sqrt(this->bins0[i]);
    }
    double scaleFactor = 0;
    for (unsigned int i = 0; i != until; ++i) {
      scaleFactor += exp(-(0.5*this->FD.externalMomentum - this->mu)*(i + 0.5)*this->dt)
                  /(sqrt(this->bins0[i]) * Z);
    }

    // greens function
    for (unsigned int i = 0; i != this->maxLength/this->dt; ++i) {
      myfile << fixed << setprecision(7) << abs(this->bins[i]*scaleFactor + singularityFix[i]);
      if (i < this->maxLength/this->dt - 1) {
        myfile << " ";
      }
    }
  }

  myfile.close();
}













double DiagrammaticMonteCarlo::Udouble (double from, double to) {
  uniform_real_distribution<double> distribution(from, to);
  return distribution(this->mt);
}

int DiagrammaticMonteCarlo::Uint (int fromIncluded, int toIncluded) {
  uniform_int_distribution<int> distribution(fromIncluded, toIncluded);
  return distribution(this->mt);
}

Vector3d DiagrammaticMonteCarlo::calculateMeanP (shared_ptr<Vertex> v1, shared_ptr<Vertex> v2) {

  Vector3d meanP(0, 0, 0);

  if (v1 != v2) {
    if (v1->G[1]->end == v2) {
      // if dt = 0
      meanP += v1->G[1]->momentum;

      if ( ! isfinite(meanP[0]) || ! isfinite(meanP[1]) || ! isfinite(meanP[2])) {
        cout << "-----------" << endl
             << "Overflow at DiagrammaticMonteCarlo::calculateMeanP" << endl
             << "meanP= " << meanP.transpose() << endl
             << "t1= " << v1->position << endl
             << "t2= " << v2->position << endl
             << "-----------" << endl;
      }
    } else {
      double dt = v2->position - v1->position;

      // loop through electrons between phonon propagator ends
      shared_ptr<Electron> g = v1->G[1];
      do {
        meanP += g->momentum*(g->end->position - g->start->position);
      } while (g->end != v2 && (g = g->end->G[1]));

      if (dt == 0 || ! isfinite(meanP[0]) || ! isfinite(meanP[1]) || ! isfinite(meanP[2])) {
        cout << "-----------" << endl
             << "Overflow at DiagrammaticMonteCarlo::calculateMeanP" << endl
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

Vector3d DiagrammaticMonteCarlo::calculateP0 (shared_ptr<Phonon> d) {
  Vector3d meanP = this->calculateMeanP(d->start, d->end);

  return meanP + d->momentum;
}

void DiagrammaticMonteCarlo::checkAcceptanceRatio (double a, string updateFunction) {
  if (abs(a - 1) > pow(10.0, -7)) {
    cout << "Acceptance ratio at " << updateFunction << ": " << setprecision(17) << a - 1 << endl;
  }
}