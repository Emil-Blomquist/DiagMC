#include "DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::parseConfig (string pathToConfigFileFromRoot) {

  string path = this->argv[0]; // program full path + name of binary file
  path.erase(path.find_last_of('/') + 1); // remove name of binary file

  string line;
  ifstream configFile(path + "../" + pathToConfigFileFromRoot);

  if (configFile.is_open()) {
    while (getline(configFile, line)) {
      istringstream currentLine(line);
      string key;
      if (getline(currentLine, key, '=')) {
        string val;
        if (getline(currentLine, val)) {


          // will print overflow exceptions etc. (code is half as fast with this activated)
          if (key == "debug_mode") {
            this->debug = stol(val);
          }
          // will print acceptance ratios etc. (debug must be true)
          if (key == "loud_mode") {
            this->loud = stol(val);
          }


          // if the diagram should have external legs or not
          if (key == "external_legs") {
            this->externalLegs = stol(val);
          }
          // if we should only take into account irreducible diagrams
          if (key == "reducible_diagrams") {
            this->reducibleDiagrams = stol(val);
          }
          // if we should sample only skeleton diagrams (reducibleDiagrams must be false)
          if (key == "skeleton_diagrams") {
            this->skeletonDiagrams = stol(val);
          }
          // if we should let the external momentum vary or not
          if (key == "fixed_external_momentum") {
            this->fixedExternalMomentum = stol(val);
          }
          // if we want to use Dyson equation (fixedExternalMomentum must be false (REALLY!?))
          // if we want to output using dyson or not
          if (key == "dyson") {
            this->Dyson = stol(val);
          }
          // wether or not we should employ boldification (fixedExternalMomentum must be false and Dyson must me set true)
          if (key == "bold") {
            this->bold = stol(val);
          }


          // for when to bin the diagram
          if (key == "minimum_diagram_order") {
            this->minDiagramOrder = stoul(val);
          }
          // raise diagrm order will look at this (zero -> turned off)
          if (key == "maximum_diagram_order") {
            this->maxDiagramOrder = stoul(val);
          }


          // if bold, how many iterations in the bold scheme shall be done
          if (key == "number_of_bold_iterations") {
            this->numBoldIterations = stoul(val);
          }
          if (key == "time_doing_mc") {
            this->numSecsDoingMCperCorePerBoldItr = stod(val);
          }
          if (key == "time_doing_dmc") {
            this->numSecsDoingDMCperCorePerBoldItr = stod(val);
          }


          // how often in seconds we should save by writing to file during each bold iteration
          if (key == "number_of_temporary_saves") {
            this->numTempDMCsaves = stoul(val);
          }


          if (key == "external_momentum") {
            this->initialExternalMomentum = Vector3d{0, 0, stod(val)};
          }
          if (key == "chemical_potential") {
            this->mu = stod(val);
          }
          if (key == "coupling_constant") {
            this->alpha = stod(val);
          }


          if (key == "dt") {
            this->dt = stod(val);
          }
          if (key == "dp") {
            this->dp = stod(val);
          }
          if (key == "max_t") {
            this->maxLength = stod(val);
          }
          if (key == "max_p") {
            this->maxMomenta = stod(val);
          }

        }
      }
    }
  }


}