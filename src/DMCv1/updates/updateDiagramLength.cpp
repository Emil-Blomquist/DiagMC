// #include "../DiagrammaticMonteCarloV1.h"

// double DiagrammaticMonteCarloV1::updateDiagramLength (double param) {

//   // old configuration value


//   double
//     l = 0.5*this->FD.end->G[0]->momentum.squaredNorm() - this->mu,
//     r = this->Udouble(0, 1),
//     tmin = this->FD.end->G[0]->start->position,
//     dt = -log(1 - r + r*exp(-l*(this->maxLength - tmin)))/l;

//   double oldVal, oldwInvt, wInvt;
//   if (this->debug) {
//     oldVal = this->FD();
//     oldwInvt = exp(l*(this->FD.end->position - tmin))*(1 - exp(-l*this->maxLength - tmin))/l;
//     wInvt = exp(l*dt)*(1 - exp(-l*this->maxLength - tmin))/l;
//   }

//   this->FD.setLength(tmin + dt);
//   this->FD.end->setPosition(tmin + dt);

//   if (this->debug) {
//     double
//       val = this->FD(),
//       a = val/oldVal * wInvt/oldwInvt;

//     if (this->loud) {
//       cout << "updateDiagramLength: " << a << endl;
//     }
//   }

//   return 1;
// }