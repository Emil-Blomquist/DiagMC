#include "DiagrammaticMonteCarlo.h"

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
         << ( this->bold ? " bolditr=" + to_string(this->boldIteration) : "" )
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

  if ( ! this->fixedExternalMomentum && this->N0) {
    // open file
    myfile.open(path + "../data/" + fileName.substr(1) + ".txt", ios_base::app);

    // the quantity of interest to us
    Array<double, Dynamic, Dynamic> quantitiyOfInterest = ArrayXXd::Zero(this->hist.rows(), this->hist.cols());

    if (this->Dyson) {
      // will become dG
      this->doDyson(quantitiyOfInterest);
      if (this->bold) {
        // save for later use
        this->dG = quantitiyOfInterest;
      }
    } else {
      // will become G
      this->normalizedHistogram(quantitiyOfInterest);
    }

    // histogram corresponding to higher order diagrams
    for (unsigned int i = 0; i != quantitiyOfInterest.rows(); i++) {
      for (unsigned int j = 0; j != quantitiyOfInterest.cols(); j++) {
        if (this->Dyson) {
          double
            p = (i + 0.5)*this->dp,
            t = (j + 0.5)*this->dt;

          myfile << fixed << setprecision(6) << quantitiyOfInterest(i, j) + exp((this->mu - 0.5*p*p)*t);
        } else {
          myfile << fixed << setprecision(6) << quantitiyOfInterest(i, j);
        }

        if (j + 1 < quantitiyOfInterest.cols()) {
          myfile << " ";
        }
      }
      if (i + 1 < quantitiyOfInterest.rows()) {
        myfile << endl;
      }
    }
  } else if (this->fixedExternalMomentum && this->bins0[0] > 0) {
    // to remove the error due the combination of a singular diagram and a discretized time
    vector<double> singularityFix((int) round(this->maxLength/this->dt), 0);
    if ( ! this->externalLegs && this->minDiagramOrder <= 1) {
      for (unsigned int i = 0; i != singularityFix.size(); ++i) {
        singularityFix[i] = this->alpha*(
          exp(-0.5*(2*i + 1)*this->dt)/sqrt(M_PI*0.5*(2*i + 1)*this->dt)
          - (erf(sqrt((i + 1)*this->dt)) - erf(sqrt(i*this->dt)))/this->dt
        );
      }
    }

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
      scaleFactor += exp(-(0.5*pow(this->FD.externalMomentum, 2) - this->mu)*(i + 0.5)*this->dt)
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