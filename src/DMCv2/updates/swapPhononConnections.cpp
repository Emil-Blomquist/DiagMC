#include "../DiagrammaticMonteCarloV2.h"

void DiagrammaticMonteCarloV2::swapPhononConnections (double param) {
  // requirement to swap: must be at least of order 2
  if (this->FD.Ds.size() < 2) {
    return;
  }

  double oldVal;
  if (this->debug) {
    oldVal = this->FD();
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
    return;
  }

  double
    c1 = (v1->D[1]) ? 1 : -1,
    c2 = (v2->D[1]) ? 1 : -1,
    t = v2->position - v1->position,
    Eafter = 0.5*(g->momentum + c1*d1->momentum - c2*d2->momentum).squaredNorm(),
    Ebefore = 0.5*g->momentum.squaredNorm();

  double exponent = -t*(Eafter - Ebefore + c2 - c1);

  // acceptance ration
  double a;
  if (exponent > 700) {
    a = 1;
  } else {
    a = exp(exponent);
  }

  // accept or reject update
  bool accepted = false;
  if (a > this->Udouble(0, 1)) {
    this->FD.swapPhonons(v1, v2);

    accepted = true;
  }

  if ( ! isfinite(a)) {
    cout
      << "-------------------------" << endl
      << "DMC::swapPhononConnections: nan encountered" << endl
      << "a=" << a << endl
      << "exponent=" << exponent << endl
      << "Eafter=" << Eafter << endl
      << "Ebefore=" << Ebefore << endl
      << "-------------------------" << endl;
  }

  if (this->debug) {
    double val = this->FD();

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::swapPhononConnections " << endl
           << "accepted=" << accepted << endl
           << "a=" << a << endl
           << "val/oldVal=" << val/oldVal << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "exponent=" << exponent << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "swapPhononConnections: "  << this->FD.Ds.size() << " " << a << endl;
    }
  }
}