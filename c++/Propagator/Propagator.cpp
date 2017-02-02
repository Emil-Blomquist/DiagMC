#include "Propagator.h"

Propagator::Propagator (int _type, Vector3d _momentum) {
  type = _type;
  momentum = _momentum;

  start = NULL;
  end = NULL;
}

void Propagator::print () {
  cout << "Propagator:" << endl
       << "\ttype: " << type << endl
       << "\tmomentum: " << momentum.transpose() << endl
       << "\tstart: " << start << endl
       << "\tend: " << end << endl;
}

void Propagator::setMomentum (Vector3d p) {
  momentum = p;
}

void Propagator::setStart (Vertex *v) {
  // unlink propagator from previous vertex
  if (start) {
    if (type == 0) {
      start->setOutgoingG(NULL);
    } else {
      start->setOutgoingD(NULL);
    }
  }

  // link propagator to new vertex
  if (type == 0) {
    v->setOutgoingG(this);
  } else {
    v->setOutgoingD(this);
  }

  // link vertex to propagator
  start = v;
}

void Propagator::setEnd (Vertex *v) {
  // unlink propagator from previous vertex
  if (end) {
    if (type == 0) {
      end->setIngoingG(NULL);
    } else {
      end->setIngoingD(NULL);
    }
  }

  // link propagator to new vertex
  if (type == 0) {
    v->setIngoingG(this);
  } else {
    v->setIngoingD(this);
  }

  // link vertex to propagator
  end = v;
}