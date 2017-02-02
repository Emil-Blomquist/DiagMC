#ifndef ELECTRON_H
#define ELECTRON_H

#include "Propagator.h"

class Electron: public Propagator {

  public:
    Electron (Vector3d _p);
};

#endif