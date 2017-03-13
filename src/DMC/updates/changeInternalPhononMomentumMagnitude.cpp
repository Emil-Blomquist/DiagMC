#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarloV2::changeInternalPhononMomentumMagnitude (double param) {
  // requirement to lower: must be at least of order 1
  if (this->FD.Ds.size() == 0) {
    return;
  }


  // select internal phonon on random
  shared_ptr<Phonon> d = this->FD.Ds[this->Uint(0, this->FD.Ds.size() - 1)];


  if (d->end->position - d->start->position == 0) {
    // if time difference is zero momentum does not matter...
    cout << "DMC::changeInternalPhononMomentumMagnitude: dt = 0 -> return" << endl;
    return;
  }


  double
    t1 = d->start->position,
    t2 = d->end->position,
    std = (t2 - t1 < pow(10.0, -10.0) ? 100000 : 1/sqrt(t2 - t1)); 
  normal_distribution<double> normal(0.0, std);
  
  double
    q = abs(normal(this->mt)),
    wInvq = sqrt(0.5*M_PI*pow(std, 2.0)) * exp(0.5*pow(q/std, 2.0));

  Vector3d
    Q(q*sin(d->theta)*cos(d->phi), q*sin(d->theta)*sin(d->phi), q*cos(d->theta)),
    oldQ = d->momentum,
    dQ = Q - oldQ,
    meanP = this->calculateMeanP(d->start, d->end);

  double
    oldq = oldQ.norm(),
    oldwInvq = sqrt(0.5*M_PI*pow(std, 2.0)) * exp(0.5*pow(oldq/std, 2.0));

  double
    dq2 = dQ.squaredNorm(),
    exponent = (dQ.dot(meanP) - 0.5*dq2)*(t2 - t1);

  double a = wInvq/oldwInvq * exp(exponent);

  double oldVal;
  if (this->debug) {
    oldVal = this->FD();
  }

  // accept or reject update
  bool accepted = false;
  if (a > this->Udouble(0, 1)) {

    // set new momentum
    this->FD.setInternalPhononMomentum(d, Q);

    accepted = true;
  }



  if (this->debug) {
    double val = this->FD();

    if (a < 0 || ! isfinite(a)) {
      cout << "--------------------------------------------------------------------" << endl
           << "overflow at DMC::changeInternalPhononMomentumMagnitude " << endl
           << "a=" << a << endl
           << "accepted=" << accepted << endl
           << "order=" << this->FD.Ds.size() << endl
           << "val=" << val << endl
           << "oldVal=" << oldVal << endl
           << "Q=" << Q.transpose() << endl
           << "--------------------------------------------------------------------" << endl;
    } else if (this->loud) {
      cout << "changeInternalPhononMomentumMagnitude " << a << endl;
    }
  }
}





// -----------
// Overflow at DiagrammaticMonteCarloV2::calculateQ
// Q= inf -nan -inf
// P0=   0.628175   0.455531 -0.0329462
// Ep=   0.808818   0.586527 -0.0424205
// Eo1=          0 -0.0721364  -0.997395
// Eo2=  -0.588059   0.806711 -0.0583452
// q= inf
// theta= 1.21672
// phi= 5.11482
// -----------
// -------------------------
// DMC::changeInternalPhononMomentumMagnitude: nan encountered
// Q= inf -nan -inf
// q=inf
// erf_inv=0.533844
// param1=0
// param2=0.269286
// param3=0.549732
// dt=0                <-----
// d->theta=1.21672
// P0=  0.628175   0.455531 -0.0329462
// p0=0.776658
// -------------------------