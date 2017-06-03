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
  if (iterationNum) {
    myfile << " Ntemp=" << iterationNum;
  }
  myfile << " --------" << endl;
  myfile.close();

  if ( ! this->fixedExternalMomentum && this->N0) {
    // open file
    myfile.open(path + "../data/" + fileName.substr(1) + ".txt", ios_base::app);

    // the quantity of interest to us
    Array<double, Dynamic, Dynamic> quantitiyOfInterest;

    // if (this->Dyson) {
    //   // will become dG
    //   this->doDyson(quantitiyOfInterest);
    // } else {
    //   // will become G
    // }

    // histogram corresponding to higher order diagrams
    for (unsigned int i = 0; i != quantitiyOfInterest.rows(); i++) {
      for (unsigned int j = 0; j != quantitiyOfInterest.cols(); j++) {
        if (this->Dyson) {
          double
            p = (i + 0.5)*this->dp,
            t = (j + 0.5)*this->dt;

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
  } else if (this->fixedExternalMomentum && this->N0) {

    // open file
    myfile.open(path + "../data/" + fileName.substr(1) + ".txt", ios_base::app);

    // times
    myfile << fixed << setprecision(6)
           << "np.linspace(0," << (this->maxLength/this->dt - 1)*this->dt << ", " << this->maxLength/this->dt << ") + " << 0.5*this->dt << endl;


    // the quantity of interest to us
    Array<double, Dynamic, 1> quantitiyOfInterest;

    this->normalizedHistogram(quantitiyOfInterest);

    // add MC contribution for S1
    for (unsigned int j = 0; j != this->S1mc.cols(); j++) {
      quantitiyOfInterest(j) += this->S1mc(j);
    }

    // histogram corresponding to higher order diagrams
    for (unsigned int j = 0; j != quantitiyOfInterest.size(); j++) {
      if (this->Dyson) {
        double
          p = this->initialExternalMomentum.norm(),
          t = (j + 0.5)*this->dt;

        // myfile << fixed << setprecision(6) << quantitiyOfInterest(i, j) + exp((this->mu - 0.5*p*p)*t);
        myfile << fixed << setprecision(6) << quantitiyOfInterest(j);
      } else {
        myfile << fixed << setprecision(6) << quantitiyOfInterest(j);
      }

      if (j + 1 < quantitiyOfInterest.size()) {
        myfile << " ";
      }
    }
  }

  myfile.close();
}