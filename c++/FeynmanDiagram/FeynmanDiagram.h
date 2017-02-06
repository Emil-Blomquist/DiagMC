#ifndef FEYNMAN_DIAGRAM_H
#define FEYNMAN_DIAGRAM_H

#include <iostream>
#include <memory>
#include <vector>
#include <sstream>

#include <eigen3/Eigen/Dense>
#include <SFML/Graphics.hpp>

using namespace Eigen;
using namespace std;

#include "../Propagator/Electron.h"
#include "../Propagator/Phonon.h"
#include "../Vertex/Vertex.h"

class FeynmanDiagram {
  private:
    // // used for plotting feynman diagram
    // int windowWidth, windowHeight, horizontalMargin, verticalMargin;
    // double thickness;
    // void drawElectron (sf::RenderWindow*, double, double, double, sf::Color);
    // void drawPhonon (sf::RenderWindow*, double, double, double, sf::Color);
    // void drawVertex (sf::RenderWindow*, double, double, bool);
    // void drawText (sf::RenderWindow*, sf::Font*, double, double, double, double);
    // sf::Color colorCode (double, double, Vector3d);
    // int longestPhonon ();

  public:
    double length, couplingConstant, chemicalPotential;
    Vector3d externalMomentum;

    shared_ptr<Vertex> start, end;

    vector<shared_ptr<Electron>> Gs;
    vector<shared_ptr<Phonon>> Ds;

    // list<Vertex> Vs;
    // list<Electron> Gs;
    // list<Phonon> Ds;

    FeynmanDiagram (Vector3d, double, double, double);

    void print ();

    shared_ptr<Vertex> insertVertex (int, double);
    shared_ptr<Phonon> addInternalPhonon (shared_ptr<Vertex>, shared_ptr<Vertex>, Vector3d, double, double);

    // void plot ();
};

#endif