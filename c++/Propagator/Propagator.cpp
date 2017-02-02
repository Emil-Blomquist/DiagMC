#include "Propagator.h"
#include "../Vertex/Vertex.h"

Propagator::Propagator (int _type, Vector3d _momentum, double _theta, double _phi) {
  type = _type;
  momentum = _momentum;
  theta = _theta;
  phi = _phi;

  start = NULL;
  end = NULL;
}

void Propagator::print () {
  cout << "Propagator:" << endl
       << "\ttype: " << type << endl
       << "\tmomentum: " << momentum.transpose() << endl
       << "\ttheta: " << theta << endl
       << "\tphi: " << phi << endl
       << "\tstart position: " << start << endl
       << "\tend position: " << end << endl;
}

void Propagator::setMomentum (Vector3d p) {
  momentum = p;
}

void Propagator::setStart (Vertex * v) {
  ////
  //// Need to unlink from previous vertex
  ////

  start = v;
}

void Propagator::setEnd (Vertex * v) {
  ////
  //// Need to unlink from previous vertex
  ////

  end = v;
}