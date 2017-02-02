#include "Phonon.h"

Phonon::Phonon (Vector3d _p, double _theta, double _phi) : Propagator(1, _p) {
  theta = _theta;
  phi = _phi;
}