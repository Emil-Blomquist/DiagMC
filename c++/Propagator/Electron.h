#ifndef ELECTRO_H
#define ELECTRO_H

#include "Propagator.h"

class Electron: public Propagator {
    Vector3d savedMomentum;
  public:
    Electron (Vector3d);
};

#endif