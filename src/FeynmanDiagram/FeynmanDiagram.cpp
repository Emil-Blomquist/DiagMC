#include "FeynmanDiagram.h"

FeynmanDiagram::FeynmanDiagram (
  Vector3d externalMomentum,
  double length,
  double couplingConstant,
  double chemicalPotential
) {
  this->externalMomentum = externalMomentum;
  this->length = length;
  this->couplingConstant = couplingConstant;
  this->chemicalPotential = chemicalPotential;

  // create vertices
  this->start.reset(new Vertex(0));
  this->end.reset(new Vertex(length));
  
  // create electronic propagator
  this->Gs.emplace_back(new Electron(externalMomentum));

  // link propagator to vertices
  this->Gs[0]->setStart(this->start);
  this->Gs[0]->setEnd(this->end);
}

double FeynmanDiagram::phononEnergy (double q) {
  return 1;
}

double FeynmanDiagram::operator() () {
  double val;

  //
  // need to change when diagram length becomes variable
  //
  val = exp(this->chemicalPotential*this->length);

  for (auto g : this->Gs) {
    val *= (*g)();
  }

  for (auto d : this->Ds) {
    val *= (*d)(this->couplingConstant, this->phononEnergy(d->q));
  }

  return val;
}

void FeynmanDiagram::setLength (double length) {
  this->length = length;
}

unsigned int FeynmanDiagram::binaryElectronSearch(shared_ptr<Electron> g, unsigned int lowerBound) {
  double position = g->end->position;

  int
    index,
    upperBound = this->Gs.size() - 1;

  do {
    index = 0.5*(upperBound + lowerBound);

    if (this->Gs[index]->end->position > position) {
      upperBound = index - 1;
    } else if (this->Gs[index]->end->position < position) {
      lowerBound = index + 1;
    } else {
      break;
    }
  } while (index != upperBound);

  if (this->Gs[index] != g) {
    // multiple vertices with the same position -> do a linear search
    cout << "FeynmanDiagram::binarySearch: forced to retreat to linear search" << endl
         << position << " vs " << this->Gs[index]->end->position << endl;

    vector<shared_ptr<Electron> >::iterator it = find(this->Gs.begin(), this->Gs.end(), g);
    index = distance(this->Gs.begin(), it);
  }

  return index;
}