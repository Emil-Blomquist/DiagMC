#ifndef PLT_H
#define PLT_H

#include <iostream>
#include <sstream>

#include <eigen3/Eigen/Dense>
#include <SFML/Graphics.hpp>
#include <cmath>

using namespace Eigen;
using namespace std;

#include "../FeynmanDiagram/FeynmanDiagram.h"

class plot {
  private:
    // used for plotting feynman diagram
    FeynmanDiagram *FD;
    int windowWidth, windowHeight, horizontalMargin, verticalMargin;
    double thickness;
    void drawElectron (sf::RenderWindow*, double, double, double, sf::Color);
    void drawPhonon (sf::RenderWindow*, double, double, double, sf::Color);
    void drawVertex (sf::RenderWindow*, double, double, bool);
    void drawText (sf::RenderWindow*, sf::Font*, double, double, double, double);
    void markVertex (sf::RenderWindow*, double, double);
    sf::Color colorCode (double, double, Vector3d);
    int longestPhonon ();

  public:
    int test;
    plot (FeynmanDiagram*);
    
    void render ();
};

#endif