#include "DiagrammaticMonteCarlo.h"

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

  this->binSize = this->maxLength/this->numBins;

  // store time at which the calculation began
  time_t rawtime;
  time (&rawtime);
  this->timeinfo = localtime(&rawtime);

  this->write2file();
  this->run();
}

void DiagrammaticMonteCarloV2::run () {
  // to reach start connfiguration
  const unsigned int untilStart = 10000000;

  // to save data under the process
  const unsigned int saveAfter = 500*1000000;

<<<<<<< HEAD

  // specify the relative probability of choosing a specific update function
  multimap<unsigned int, void (DiagrammaticMonteCarloV2::*)(double)> temp = {
    {49, &DiagrammaticMonteCarloV2::shiftVertexPosition},
    {46, &DiagrammaticMonteCarloV2::swapPhononConnections},
    {26, &DiagrammaticMonteCarloV2::changeInternalPhononMomentumDirection},
    {66, &DiagrammaticMonteCarloV2::changeInternalPhononMomentumMagnitude},
    {57, &DiagrammaticMonteCarloV2::raiseOrder}, // <- These two must have the same probability
    {57, &DiagrammaticMonteCarloV2::lowerOrder}, // <-
    {11, &DiagrammaticMonteCarloV2::changeDiagramLength}
  };

  // vector which is going to contain the specified quantity of update functions
  vector<void (DiagrammaticMonteCarloV2::*)(double)> updateMethods;

  // populate vector
  for (auto updateMethod = temp.begin(); updateMethod != temp.end(); updateMethod++) {
    for (unsigned int i = 0; i != updateMethod->first; i++) {
      updateMethods.push_back(updateMethod->second);
    }
  }

  // vector of pointers to member function of Phonon
  // vector<void (DiagrammaticMonteCarloV2::*)(double)> updateMethods = {
  //   &DiagrammaticMonteCarloV2::shiftVertexPosition, // <- 2
  //   &DiagrammaticMonteCarloV2::swapPhononConnections, // <- 1
  //   &DiagrammaticMonteCarloV2::changeInternalPhononMomentumDirection, // <- 3
  //   &DiagrammaticMonteCarloV2::changeInternalPhononMomentumMagnitude, // <- 4
  //   &DiagrammaticMonteCarloV2::raiseOrder,
  //   &DiagrammaticMonteCarloV2::lowerOrder,
  //   &DiagrammaticMonteCarloV2::changeDiagramLength
  // };

  // bins and keys for counting
  this->keys = vector<double>(this->numBins, 0);
  this->bins = vector<int>(this->numBins, 0);
=======
  // bins for counting
  this->bins = vector<unsigned long int>(this->numBins, 0);
  this->bins0 = vector<unsigned long int>(this->numBins, 0);

  // specify the relative probability of choosing a specific update function
  multimap<unsigned int, void (DiagrammaticMonteCarloV2::*)(double)> updateMethods = {
    {2, &DiagrammaticMonteCarloV2::shiftVertexPosition},
    {2, &DiagrammaticMonteCarloV2::swapPhononConnections},
    {2, &DiagrammaticMonteCarloV2::changeInternalPhononMomentumDirection},
    {2, &DiagrammaticMonteCarloV2::changeInternalPhononMomentumMagnitude},
    {2, &DiagrammaticMonteCarloV2::raiseOrder}, // <- These two must have the same probability
    {2, &DiagrammaticMonteCarloV2::lowerOrder}, // <-
    {2, &DiagrammaticMonteCarloV2::changeDiagramLength},
    {1, &DiagrammaticMonteCarloV2::changeDiagramLengthComplex}
  };

  // vector which is going to contain the specified quantity of update functions
  vector<void (DiagrammaticMonteCarloV2::*)(double)> chooseUpdateMethod;
>>>>>>> Easy-mode

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

  // Display disp(&this->FD);
  // disp.render();

  // main loop
  for (unsigned long int i = 0; i < this->numIterations; ++i) {
    auto updateMethod = chooseUpdateMethod[this->Uint(0, chooseUpdateMethod.size() - 1)];
    (this->*updateMethod)(this->param);

    if (i%saveAfter == saveAfter - 1) {
      // save temporary result
      this->write2file(i + 1);
    }

    // bin
    unsigned int index = this->FD.length/this->binSize;
    bins[index]++;

    // if zeroth order, bin again
    if (this->FD.Ds.size() == 0) {
      bins0[index]++;
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