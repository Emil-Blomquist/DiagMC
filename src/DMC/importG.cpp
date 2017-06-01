#include "DiagrammaticMonteCarlo.h"

template<typename Out>
void split(const string &s, char delim, Out result) {
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
        *(result++) = item;
    }
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, back_inserter(elems));
    return elems;
}

void DiagrammaticMonteCarlo::importG (string fileName) {

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
        vector<string> args = split(line, ' ');

        // read in important variables
        for (string arg : args) {
          vector<string> keyAndValue = split(arg, '=');

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

        this->dG = Array<double, Dynamic, Dynamic>::Zero(Np, Nt);
        this->dE = Array<double, Dynamic, Dynamic>::Zero(Np, Nt);
      } else {
        double
          p = (i + 0.5)*this->dp,
          E0 = 0.5*p*p;

        vector<string> row = split(line, ' ');
        for (unsigned int j = 0; j != this->dG.cols(); j++) {
          double t = (j + 0.5)*this->dt;
          this->dG(i, j) = stod(row[j]) - exp((this->mu - E0)*t);
        }

        i++;
      }
    }

    myfile.close();

    // cout << this->dG << endl;

    this->calculateEnergyDiff();

    // cout << this->dE << endl;



  } else {
    cout << "Unable to open file" << endl; 
  }
}