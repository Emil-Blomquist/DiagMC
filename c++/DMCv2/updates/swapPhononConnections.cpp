#include "../DiagrammaticMonteCarloV2.h"

double DiagrammaticMonteCarloV2::swapPhononConnections (double param) {
  // requirement to swap: must be at least of order 2
  if (this->FD.Ds.size() < 2) {
    return 0;
  }

  // select internal electron propagator on random
  shared_ptr<Electron> g = this->FD.Gs[this->Uint(1, this->FD.Gs.size() - 2)];
  
  shared_ptr<Vertex>
    v1 = g->start,
    v2 = g->end;

  // fetch connected phonons
  shared_ptr<Phonon>
    d1 = v1->D[0] ? v1->D[0] : v1->D[1],
    d2 = v2->D[0] ? v2->D[0] : v2->D[1];

  // we do not allow for phonon to be reversed
  if (d1 == d2) {
    return 0;
  }

  double
    c1 = (v1->D[1]) ? 1 : -1,
    c2 = (v2->D[1]) ? 1 : -1,
    t = v2->position - v1->position,
    Eafter = 0.5*(g->momentum + c1*d1->momentum - c2*d2->momentum).squaredNorm(),
    Ebefore = 0.5*g->momentum.squaredNorm();

  double exponent = -t*(Eafter - Ebefore + c2 - c1);

  // to prevent overflow error
  if (exponent > 700) { exponent = 700; }

  // acceptance ratio
  double a = exp(exponent);

  double oldVal;
  if (this->debug) {
    oldVal = this->FD();
  }

  this->FD.swapPhonons(v1, v2);

  if (this->debug) {
    double val = this->FD();

    if (exponent == 700) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow prevented at DMC::swapPhononConnections " << a*oldVal/val << " val=" << val << " oldVal=" << oldVal << endl
           << "--------------------------------------------------------------------" << endl;
    } else {
      if (this->loud) { cout << "swapPhononConnections: " << a*oldVal/val << endl; }
    }
  }

  return a;
}