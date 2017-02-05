#ifndef FEYNMAN_DIAGRAM_H
#define FEYNMAN_DIAGRAM_H

#include <iostream>
#include <list>
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
    // used for plotting feynman diagram
    int windowWidth, windowHeight, horizontalMargin, verticalMargin;
    double thickness;
    void drawElectron (sf::RenderWindow*, double, double, double);
    void drawPhonon (sf::RenderWindow*, double, double, double);
    void drawVertex (sf::RenderWindow*, double, double, bool);
    void drawText (sf::RenderWindow*, sf::Font*, double, double, double, double);
    int longestPhonon ();

  public:
    double length, couplingConstant, chemicalPotential;
    Vector3d externalMomentum;

    list<Vertex> Vs;
    list<Electron> Gs;
    list<Phonon> Ds;

    FeynmanDiagram (Vector3d, double, double, double);

    void print ();

    Vertex *insertVertex (list<Vertex>::iterator, double);
    Phonon *addInternalPhonon (list<Vertex>::iterator, list<Vertex>::iterator, Vector3d, double, double);

    void plot ();
};

#endif