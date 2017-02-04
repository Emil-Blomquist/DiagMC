#ifndef PHONON_H
#define PHONON_H

#include "Propagator.h"

class Phonon: public Propagator {
    Vector3d savedMomentum;
    // double savedTheta, savedPhi;
  public:

    double theta, phi;

    Phonon (Vector3d _p, double _theta, double _phi);

    void print ();

    void setStart (Vertex *);
    void setEnd (Vertex *);
};

#endif