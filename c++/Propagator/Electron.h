#ifndef ELECTRON_H
#define ELECTRON_H

#include "Propagator.h"

class Electron: public Propagator, public enable_shared_from_this<Electron> {

  public:
    Electron (Vector3d);

    void print ();

    void setStart (shared_ptr<Vertex>);
    void setEnd (shared_ptr<Vertex>);
};

#endif