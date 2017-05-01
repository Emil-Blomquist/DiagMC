#include "FeynmanDiagram.h"

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

// void FeynmanDiagram::insertIntoHashTable (shared_ptr<Electron> g) {
//   // decrease precision so numerical errors wont affect
//   const double key = round(g->momentum[2] * this->tolerance) / this->tolerance;

//   this->electronHashTable.insert({key, g});
// }

// void FeynmanDiagram::removeFromHashTable (shared_ptr<Electron> g) {
//   // decrease precision so numerical errors wont affect
//   const double key = round(g->momentum[2] * this->tolerance) / this->tolerance;

//   // obtain iterator pointer pointing towards the propagator
//   auto range = this->electronHashTable.equal_range(key);

//   if (distance(range.first, range.second) == 1) {
//     // only one propagator with this key
//     this->electronHashTable.erase(range.first);
//   } else {
//     // multiple propagators with this key

//     // linear search
//     for (auto i = range.first; i != range.second; i++) {
//       if (i->second == g) {
//         this->electronHashTable.erase(i);
//         return;
//       }
//     }

//     cout << "PROPAGATOR NOT FOUND" << endl;
//   }
// }


bool FeynmanDiagram::diagramIsIrreducible () {
  if (this->newStructure) {
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

// void FeynmanDiagram::printHashTable () {


//   if (this->Gs.size() - this->electronHashTable.size() != 0) {
//     cout << "Mismatch: " << this->Gs.size() << " vs " << this->electronHashTable.size() << endl;
//   }

//   // vector<double> uniqueKeys;

//   // cout << "--------------------- " << this->Gs.size() << endl;
//   // for (auto& x : this->electronHashTable) {
//   //   if (! uniqueKeys.size() || find(uniqueKeys.begin(), uniqueKeys.end(), x.first) == uniqueKeys.end()) {
//   //     uniqueKeys.push_back(x.first);

//   //     cout.precision(17);
//   //     cout << "\t" << fixed << x.first << ": " << this->electronHashTable.count(x.first) << endl;
//   //   }

//   // }
//   // cout << "---------------------" << endl;
// }

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