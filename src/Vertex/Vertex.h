#ifndef VERTEX_H
#define VERTEX_H

#include <iostream>
#include <memory>
#include <math.h>

using namespace std;

class Electron;
class Phonon;

class Vertex : public enable_shared_from_this<Vertex> {
  public:
    double position;

    shared_ptr<Electron> G[2];
    shared_ptr<Phonon> D[2];

    Vertex (double);

    void
      setG (int, shared_ptr<Electron>),
      setD (int, shared_ptr<Phonon>),
      setPosition (double);
};

#endif