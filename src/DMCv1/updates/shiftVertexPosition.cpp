// #include "../DiagrammaticMonteCarloV1.h"

// double DiagrammaticMonteCarloV1::shiftVertexPosition (double param) {
//   // requirement to lower: must be at least of order 1
//   if (this->FD.Ds.size() == 0) {
//     return 0;
//   }

//   // select vertex on random
//   shared_ptr<Electron> g = this->FD.Gs[this->Uint(0, this->FD.Gs.size() - 2)];
//   shared_ptr<Vertex> v = g->end;

//   // fetch available time interval
//   double
//     t1 = g->start->position,
//     t2 = v->G[1]->end->position;

//   // calculate exponent
//   double
//     c = (v->D[0]) ? -1 : 1,
//     dE = 0.5*v->G[0]->momentum.squaredNorm() - 0.5*v->G[1]->momentum.squaredNorm() - c,
//     dtdE = (t2 - t1)*dE;

//   // sample new t
//   double t, r = this->Udouble(0, 1);
//   if (-dtdE > 100) {
//     // to avoid overflow due to exponential
//     t = t2 - log(r)/dE;
//   } else {
//     t = t1 - log(1 - r*(1 - exp(-dtdE)))/dE;
//   }

//   double tOld, oldVal;
//   if (this->debug) {
//     tOld = v->position,
//     oldVal = this->FD();
//   }

//   this->FD.setVertexPosition(v, t);

//   if (this->debug) {
//     double 
//       val = this->FD(),
//       a = exp(log(val) - log(oldVal) + dE*(t - tOld));

//     if (this->loud) { cout << "shiftVertexPosition: " << a << endl; }
//   }

//   return 1;
// }