#include "MonteCarlo.h"

template<typename Out>
void split_MC(const string &s, char delim, Out result) {
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
        *(result++) = item;
    }
}

vector<string> split_MC(const string &s, char delim) {
    vector<string> elems;
    split_MC(s, delim, back_inserter(elems));
    return elems;
}

void MonteCarlo::importG (string fileName) {

  // obtain path relative binary file
  string path = this->argv[0]; // program full path + name of binary file
  path.erase(path.find_last_of('/') + 1); // remove name of binary file



  bool firstLine = true;
  unsigned int i = 0;

  string line;
  ifstream myfile(path + "../Gs/" + fileName);
  if (myfile.is_open()) {
    while (getline(myfile, line)) {

      if (firstLine) {
        firstLine = false;

        // split first line on white spaces
        vector<string> args = split_MC(line, ' ');

        // read in important variables
        for (string arg : args) {
          vector<string> keyAndValue = split_MC(arg, '=');

          // overwrite parameters
          if (keyAndValue[0] == "a") {
            this->alpha = stod(keyAndValue[1]);
          } else if (keyAndValue[0] == "mu") {
            this->mu = stod(keyAndValue[1]);
          } else if (keyAndValue[0] == "tmax") {
            this->maxLength = stod(keyAndValue[1]);
          } else if (keyAndValue[0] == "dt") {
            this->dt = stod(keyAndValue[1]);
          } else if (keyAndValue[0] == "pmax") {
            this->maxMomenta = stod(keyAndValue[1]);
          } else if (keyAndValue[0] == "dp") {
            this->dp = stod(keyAndValue[1]);
          }
        }

        // initiate dE and dG
        unsigned int
          Np = this->maxMomenta/this->dp,
          Nt = this->maxLength/this->dt;

        this->Gfull = Array<double, Dynamic, Dynamic>::Zero(Np, Nt);
        this->dG = Array<double, Dynamic, Dynamic>::Zero(Np, Nt);
        this->dE = Array<double, Dynamic, Dynamic>::Zero(Np, Nt);
      } else {
        double
          p = (i + 0.5)*this->dp,
          E0 = 0.5*p*p;

        vector<string> row = split_MC(line, ' ');
        for (unsigned int j = 0; j != this->dG.cols(); j++) {
          double t = (j + 0.5)*this->dt;

          this->Gfull(i, j) = stod(row[j]);
          this->dG(i, j) = stod(row[j]) - exp((this->mu - E0)*t);
        }

        i++;
      }
    }

    myfile.close();

    // calcualte dE
    this->calculateEnergyDiff();


      // if (false) {
      //   cout << "dE" << " = np.array([";
      //   for (unsigned int i = 0; i != this->dE.cols(); i++) {
      //     cout << this->dE(0, i);

      //     if (i < this->dE.cols() - 1) {
      //       cout << ", ";
      //     }
      //   }
      //   cout << "])" << endl << endl;


      //   cout << "dG" << " = np.array([";
      //   for (unsigned int i = 0; i != this->dG.cols(); i++) {
      //     cout << this->dG(0, i);

      //     if (i < this->dG.cols() - 1) {
      //       cout << ", ";
      //     }
      //   }
      //   cout << "])" << endl << endl;
      // } else {
      //   unsigned int it = round(1/this->dt);

      //   cout << it << endl;

      //   cout << "dE" << " = np.array([";
      //   for (unsigned int i = 0; i != this->dE.rows(); i++) {
      //     cout << this->dE(i, it);

      //     if (i < this->dE.rows() - 1) {
      //       cout << ", ";
      //     }
      //   }
      //   cout << "])" << endl << endl;


      //   cout << "dG" << " = np.array([";
      //   for (unsigned int i = 0; i != this->dG.rows(); i++) {
      //     cout << this->dG(i, it);

      //     if (i < this->dG.rows() - 1) {
      //       cout << ", ";
      //     }
      //   }
      //   cout << "])" << endl << endl;
      // }







  } else {
    cout << "Unable to open file" << endl; 
  }
}

void MonteCarlo::calculateEnergyDiff () {
  for (unsigned int i = 0; i != this->dE.rows(); i++) {
    double
      p = (i + 0.5)*this->dp,
      E0 = 0.5*p*p;

    for (unsigned int j = 0; j != this->dE.cols(); j++) {
      double
        t = (j + 0.5)*this->dt,
        G0 = exp((this->mu - E0)*t),
        logOfG = (this->mu - E0)*t + (isfinite(this->dG(i, j)/G0) ? log(1 + abs(this->dG(i, j)/G0)) : 0);

      this->dE(i, j) = this->mu - E0 - logOfG/t;

      if ( ! isfinite(this->dE(i, j))) {
        cout << p << " " << t << fixed << setprecision(16) << ": " << log(G0 + this->dG) << " vs "
             << (this->mu - 0.5*p*p)*t + log(1 + this->dG/G0) << " but " << this->dG/G0 << endl;
      }
    }
  }
}