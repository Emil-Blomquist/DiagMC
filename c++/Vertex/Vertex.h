#ifndef VERTEX_H
#define VERTEX_H

#include <iostream>

using namespace std;

class Electron;
class Phonon;

class Vertex {
  public:
    double position;

    Electron *G[2];
    Phonon *D[2];

    Vertex (double);

    void print ();

    void setIngoingG (Electron *prop) {this->G[0] = prop;}
    void setOutgoingG (Electron *prop) {this->G[1] = prop;}
    void setIngoingD (Phonon *prop) {this->D[0] = prop;}
    void setOutgoingD (Phonon *prop) {this->D[1] = prop;}

    // void setPosition (double);
};

#endif