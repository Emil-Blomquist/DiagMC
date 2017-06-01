#include "MonteCarlo.h"
MonteCarlo::MonteCarlo (
  Vector3d externalMomentum, 
  double alpha,
  double mu,
  unsigned long int numIterations,
  double tStart,
  double tEnd,
  char **argv
) {
  // seed random generator
  unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
  this->mt.seed(seed1);
  
  // store current directory
  this->argv = argv;

  // store time at which the calculation began
  time_t rawtime;
  time (&rawtime);
  this->timeinfo = localtime(&rawtime);

  // obtain path relative binary file
  // program full path + name of binary file
  this->path = argv[0];
  // remove name of binary file
  this->path.erase(this->path.find_last_of('/') + 1);

  this->externalMomentum = externalMomentum;
  this->mu = mu;
  this->alpha = alpha;
  this->tStart = tStart;
  this->tEnd = tEnd;
  this->numIterations = numIterations;


  this->dtOut = 0.02;
  this->externalLegs = false;
  this->irreducibleDiagrams = true;
  this->maxOrder = 1;

  // Doing bold scheme
  if (true) {
    this->importG("Glarge.txt");
  }

  // create times vector
  for (unsigned int i = tStart/dtOut; i <= tEnd/dtOut - 1; i++) {
    this->times.push_back((i + 0.5)*dtOut);
  }

  // initiate container storing the values of the Greens function
  this->values = {};

  this->write2file();
  this->run();
}

void MonteCarlo::run () {
  for (unsigned int i = 0; i != this->times.size(); i++) {
    // initiate
    this->values.push_back(0);

    // zeroth order
    if ( ! this->irreducibleDiagrams) {
      this->values[i] = this->G(this->externalMomentum, 0, this->times[i]);
    }

    if (this->maxOrder >= 1) this->diagramOrder1(this->times[i], i);

    if (this->maxOrder >= 2) this->diagramOrder2(this->times[i], i);

    this->write2file();
  }
}



void MonteCarlo::write2file () {
  // create date and time string
  char buffer[80];
  strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", this->timeinfo);
  string dateAndTimeString(buffer);

  // create file name
  stringstream stream;
  stream << fixed << setprecision(7)
         << "p=" << this->externalMomentum.norm()
         << " a=" << this->alpha
         << " mu=" << this->mu
         << " tstart=" << this->tStart
         << " tend=" << this->tEnd
         << " dt=" << this->dtOut
         << " N=" << this->numIterations
         << " date=" << dateAndTimeString;
  string fileName = stream.str();

  // create data folder if not already existing
  system(("mkdir -p " + this->path + "../data").c_str());

  // write to file
  ofstream myfile;
  myfile.open(this->path + "../data/" + fileName + ".txt");
  myfile << "-------- " + fileName + " --------" << endl;

  // write data
  for (unsigned int i = 0; i != this->values.size(); i++) {
    myfile << fixed << setprecision(7) << this->times[i];
    if (i < this->values.size() - 1) {
      myfile << " ";
    }
  }
  myfile << endl;
  for (unsigned int i = 0; i != this->values.size(); i++) {
    myfile << fixed << setprecision(7) << this->values[i];
    if (i < this->values.size() - 1) {
      myfile << " ";
    }
  }

  myfile.close();
}


double MonteCarlo::G (Vector3d P, double t1, double t2) {
  double
    dG = 0,
    dE = 0,
    p = P.norm(),
    E = 0.5*p*p;

  if (this->dG.size()) {
    unsigned int
      pi = min(p/this->dp, this->dG.rows() - 1.0),
      ti = (t2 - t1)/this->dt;

    dG = this->dG(pi, ti);
    dE = this->dE(pi, ti);
  }

  return exp((this->mu - E - dE)*(t2 - t1));
  // return exp((this->mu - E)*(t2 - t1)) + dG;
}

double MonteCarlo::D (Vector3d Q, double theta, double t1, double t2) {
  double E = this->phononDispersionRelation(Q);

  return this->alpha/(sqrt(8)*M_PI*M_PI) * sin(theta) * exp(-E*(t2 - t1));
}

double MonteCarlo::phononDispersionRelation (Vector3d Q) {
  return 1;
}


double MonteCarlo::Udouble (double from, double to) {
  uniform_real_distribution<double> distribution(from, to);
  return distribution(this->mt);
}

double MonteCarlo::Ndouble (double std) {
  normal_distribution<double> distribution(0.0, std);
  return distribution(this->mt);
}

