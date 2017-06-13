#include "FeynmanDiagram.h"

FeynmanDiagram::FeynmanDiagram () {
  // default constructor
}

FeynmanDiagram::FeynmanDiagram (
  Vector3d ExternalMomentum,
  double length,
  double couplingConstant,
  double chemicalPotential
) {
  this->ExternalMomentum = ExternalMomentum;
  this->externalMomentum = ExternalMomentum.norm();
  this->length = length;
  this->couplingConstant = couplingConstant;
  this->chemicalPotential = chemicalPotential;

  // initial value
  this->isIrreducible = false;
  this->newStructure = true;

  // tolerance for hash table
  // this->tolerance = 10000000;

  // create vertices
  this->start.reset(new Vertex(0));
  this->end.reset(new Vertex(length));
  
  // create electronic propagator
  this->Gs.emplace_back(new Electron(ExternalMomentum));

  // add the default electron propagator to the hash table
  // this->insertIntoHashTable(this->Gs[0]);

  // link propagator to vertices
  this->Gs[0]->setStart(this->start);
  this->Gs[0]->setEnd(this->end);
}

double FeynmanDiagram::phononEnergy (double q) {
  // OBS:  1) Not implemented in update functions 
  //       2) need to modify MC/DMC boundary if this becomes non-constant
  return 1;
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



bool FeynmanDiagram::isIrreducibleDiagram (bool checkAnyway) {
  if (this->newStructure || checkAnyway) {
    // need to check since structure is changed
    map<shared_ptr<Phonon>, bool> phononMap;

    shared_ptr<Vertex> v = this->start;
    do {
      if (v->D[1]) {
        // outgoing phonon
        phononMap[v->D[1]] = true;
      } else if (v->D[0]) {
        // ingoing phonon
        phononMap.erase(v->D[0]);
        if (phononMap.empty()) {
          if (v != this->end) {
            this->isIrreducible = false;
          } else {
            this->isIrreducible = true;
          }
          break;
        }
      }
    } while (v != this->end && (v = v->G[1]->end));

    this->newStructure = false;
  }
  
  return this->isIrreducible;
}

bool FeynmanDiagram::phononGoesBeyond (shared_ptr<Phonon> d, shared_ptr<Vertex> end) {
  // compare times
  if (d->end->position < end->position) {
    return false;
  } else if (d->end->position > end->position) {
    return true;
  } else {
    // times are equal, we need to do a linear search

    shared_ptr<Vertex> v = end;
    do {
      if (v == d->end) {
        // if the end is found between "d->start" and the end of the subdiagram
        return false;
      }
    } while ((v = v->G[0]->start) && v != d->start);

    return true;
  }
}

bool FeynmanDiagram::isIrreducibleSubdiagram (shared_ptr<Vertex> start, shared_ptr<Vertex> end) {
  // start and end of the subdiagram

  // need to check since structure is changed
  set<shared_ptr<Phonon>> phononSet;

  shared_ptr<Vertex> v = start;
  do {
    if (v->D[1]) {
      // outgoing phonon

      // check if this phonon preceedes the end vertex
      // if this is the case, clear the hash table since they intersect each other
      // otherwise add to the hashtable
      if (this->phononGoesBeyond(v->D[1], end)) {
        phononSet.clear();
      } else {
        phononSet.insert(v->D[1]);
      }
    } else if (v->D[0]) {
      // ingoing phonon

      // check if this phonon exists in the hash table
      // if this is the case, simply remove and check wheter or not the hash table is empty
      // if not, clear the hash table since they intersect each other
      if (phononSet.find(v->D[0]) != phononSet.end()) {
        phononSet.erase(v->D[0]);
        if (phononSet.empty()) {
          return false;
        }
      } else {
        phononSet.clear();
      }

    }
  } while (v != end && (v = v->G[1]->end));
  
  // if we have reached this line the subdiagram must be irreducible
  return true;
}

bool FeynmanDiagram::isSkeletonDiagram () {
  if (this->newStructure) {
    this->newStructure = false;

    // the first order diagram is an exception
    if (this->Ds.size() == 1) {
      this->isSkeleton = true;
      return true;
    }

    // first we should make sure that the diagram is not reducible
    if ( ! this->isIrreducibleDiagram(true)) {
      this->isSkeleton = false;
      return false;
    }

    // next we need to make sure that each subdiagram also is irreducible
    for (shared_ptr<Phonon> d : this->Ds) {
      if (
        d->start->G[1]->end == d->end ||
        ! this->isIrreducibleSubdiagram(d->start->G[1]->end, d->end->G[0]->start)
      ) {
        this->isSkeleton = false;
        return false;
      }
    }

    // if we reach this line, the diagram sure is a skeleton one
    this->isSkeleton = true;
  }

  return this->isSkeleton;
}

string FeynmanDiagram::diagramName () {
  string name = "";

  map<shared_ptr<Phonon>, unsigned int> hash;

  unsigned int n = 1;

  shared_ptr<Vertex> v = this->start;
  do {
    if (v->D[1]) {
      // outgoing phonon
      hash[v->D[1]] = n;
      name += to_string(n++);
    } else if (v->D[0]) {
      // ingoing phonon
      name += to_string(hash.find(v->D[0])->second);
    }
  } while (v != this->end && (v = v->G[1]->end));


  if (name == "") {
    return "0";
  }

  return name;
}







void FeynmanDiagram::setExternalMomentum (Vector3d P, double p, Vector3d dP) {
  this->ExternalMomentum = P;
  this->externalMomentum = p;

  shared_ptr<Electron> g = this->Gs[0];
  do {
    // remove from hash table
    // this->removeFromHashTable(g);  

    g->addMomentum(dP);

    // reinsert into hash table
    // this->insertIntoHashTable(g);
  } while (g->end != this->end && (g = g->end->G[1]));
}

void FeynmanDiagram::setNewStructure () {
  this->newStructure = true;
}