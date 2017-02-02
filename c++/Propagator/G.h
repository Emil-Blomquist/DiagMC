#ifndef G_H
#define G_H

#include "Propagator.h"

class G: public Propagator {

  public:
    G (Vector3d _p, double _theta, double _phi);
};

#endif