#ifndef VERTEX_H
#define VERTEX_H

#include <iostream>
#include <memory>

using namespace std;

class Electron;
class Phonon;

class Vertex : public enable_shared_from_this<Vertex> {
    double savedPosition;
    shared_ptr<Electron> savedG[2];
    shared_ptr<Phonon> savedD[2];
  public:
    double position;

    shared_ptr<Electron> G[2];
    shared_ptr<Phonon> D[2];

    Vertex (double);

    void print ();

    void setG (int, shared_ptr<Electron>);
    void setD (int, shared_ptr<Phonon>);

    void setPosition (double);

    void save ();
    void revert ();
};

#endif