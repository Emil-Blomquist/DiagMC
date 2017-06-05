#include "DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::write2file (
  Array<unsigned long long int, Dynamic, Dynamic>& hist,
  unsigned long long int& N0,
  unsigned long long int& currNumItr
) {
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
         << " N=" << this->numIterations
         << ( this->bold ? " bolditr=" + to_string(this->boldIteration) : "" )
         << ( this->bold ? " maxn=" + to_string(this->maxDiagramOrder) : "" )
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
  if (this->numIterations != currNumItr) {
    myfile << " Ntemp=" << currNumItr;
  }
  myfile << " --------" << endl;

  if (this->N0) {
    // append the data

    // the quantity of interest to us
    Array<double, Dynamic, Dynamic> quantitiyOfInterest;
    this->normalizeHistogram(hist, N0, quantitiyOfInterest);

    // if (this->Dyson) {
    //   // will become dG
    //   this->doDyson(quantitiyOfInterest);
    // } else {
    //   // will become G
    // }

    if (this->fixedExternalMomentum) {
      // write times
      myfile << fixed << setprecision(6)
             << "np.linspace(0,"
             << (this->maxLength/this->dt - 1)*this->dt
             << ", " << this->maxLength/this->dt << ") + "
             << 0.5*this->dt << endl;
    }

    // write histogram
    for (unsigned int i = 0; i != quantitiyOfInterest.rows(); i++) {
      // double p = (this->fixedExternalMomentum ? initialExternalMomentum.norm() : (i + 0.5)*this->dp);
      for (unsigned int j = 0; j != quantitiyOfInterest.cols(); j++) {
        // double t = (j + 0.5)*this->dt;
        if (this->Dyson) {
          // myfile << fixed << setprecision(6) << quantitiyOfInterest(i, j) + exp((this->mu - 0.5*p*p)*t);
          myfile << fixed << setprecision(6) << quantitiyOfInterest(i, j);
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
  }

  myfile.close();
}