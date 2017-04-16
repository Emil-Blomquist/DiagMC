#include "DiagrammaticMonteCarlo.h"

DiagrammaticMonteCarlo::DiagrammaticMonteCarlo (
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
  unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
  this->mt.seed(seed1);

  // will print overflow exceptions etc. (code is half as fast with this activated)
  this->debug = false;
  // will print acceptance ratios etc. (debug must be true)
  this->loud = false;

  // if the diagram should have external legs or not
  this->externalLegs = false;
  this->irreducibleDiagrams = true;

  this->mu = mu;
  this->alpha = alpha;
  this->maxLength = maxLength;
  this->numIterations = numIterations;
  this->numBins = numBins;
  this->param = param;
  this->argv = argv;

  this->binSize = this->maxLength/this->numBins;

  // store time at which the calculation began
  time_t rawtime;
  time (&rawtime);
  this->timeinfo = localtime(&rawtime);

  this->write2file();
  this->run();
}

void DiagrammaticMonteCarlo::run () {
  // to reach start connfiguration
  const unsigned int untilStart = 10000000;

  // to save data under the process
  const unsigned int saveAfter = 500*1000000;

  // bins for counting
  this->bins = vector<unsigned long int>(this->numBins, 0);
  this->bins0 = vector<unsigned long int>(this->numBins, 0);

  // specify the relative probability of choosing a specific update function
  multimap<unsigned int, void (DiagrammaticMonteCarlo::*)(double)> updateMethods = {
    {10, &DiagrammaticMonteCarlo::shiftVertexPosition},
    {3, &DiagrammaticMonteCarlo::swapPhononConnections},
    {5, &DiagrammaticMonteCarlo::changeInternalPhononMomentumDirection},
    {7, &DiagrammaticMonteCarlo::changeInternalPhononMomentumMagnitude},
    {2, &DiagrammaticMonteCarlo::raiseOrder}, // <- These two must have the same probability
    {2, &DiagrammaticMonteCarlo::lowerOrder}, // <-
    {5, &DiagrammaticMonteCarlo::changeDiagramLength},
    {1, &DiagrammaticMonteCarlo::changeDiagramLengthComplex}
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

    if (
      this->externalLegs ||
      ( ! this->externalLegs && this->FD.Ds.size() >= 2)
    ) {
      if (
        ! this->irreducibleDiagrams ||
        (this->irreducibleDiagrams && this->FD.diagramIsIrreducible(this->externalLegs))
      ) {
        unsigned int index = this->FD.length/this->binSize;
        bins[index]++;
      }
    }

    // bin
    // if (this->FD.Ds.size() >= 2 && this->FD.diagramIsIrreducible(this->externalLegs)) {


        
    //   // to verify
    //   bool isIrreducibleVerified = true;
    //   map<shared_ptr<Phonon>, bool> hash;
    //   shared_ptr<Vertex> v = this->FD.start;
    //   do {
    //     if (v->D[1]) {
    //       // outgoing phonon
    //       hash[v->D[1]] = true;
    //     } else if (v->D[0]) {
    //       // ingoing phonon
    //       hash.erase(v->D[0]);
    //       if (hash.size() == 0 && v != this->FD.end) {
    //         isIrreducibleVerified = false;
    //         break;
    //       }
    //     }
    //   } while (v != this->FD.end && (v = v->G[1]->end));

    //   if ( ! isIrreducibleVerified) {
    //     cout.precision(17);

    //     cout << "actual: " << round(this->FD.externalMomentum[2] * 1000000000) / 1000000000 << endl;

    //     for (auto iter = this->FD.electronHashTable.begin(); iter != this->FD.electronHashTable.end(); ++iter) {
    //       cout << "keys: " << iter->first << endl;
    //     }

    //     Display disp(&this->FD);
    //     disp.render();
    //   }
    // }

    // if zeroth order, bin again
    if (this->FD.Ds.size() == 0) {
      unsigned int index = this->FD.length/this->binSize;
      bins0[index]++;
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
  stream << fixed << setprecision(7)
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

  if (! this->bins.empty()) {
    // data to write to file
    myfile.open(path + "../data/" + fileName + ".txt", ios_base::app);

    // times
    for (unsigned int i = 0; i != this->numBins; ++i) {
      myfile << fixed << setprecision(7) << (i + 0.5)*this->binSize;
      if (i == this->numBins - 1) {
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
      scaleFactor += exp(-(0.5*this->FD.externalMomentum.squaredNorm() - this->mu)*(i + 0.5)*this->binSize)
                  /(sqrt(this->bins0[i]) * Z);
    }

    // greens function
    for (unsigned int i = 0; i != this->numBins; ++i) {
      myfile << fixed << setprecision(7) << (double) this->bins[i]*scaleFactor;
      if (i < this->numBins - 1) {
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