#ifndef PHONON_H
#define PHONON_H

#include "Propagator.h"

class Phonon: public Propagator, public enable_shared_from_this<Phonon> {
    double savedTheta, savedPhi;
  public:
    double theta, phi;

    Phonon (Vector3d _p, double _theta, double _phi);

    void print ();

    void setStart (shared_ptr<Vertex>);
    void setEnd (shared_ptr<Vertex>);
    void setTheta (double);
    void setPhi (double);

    double operator() (double);

    void save ();
    void revert ();
};

#endif