#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarlo::changeDiagramOrder (double param) {
  
  // to make alpha independent and even out the n dependence
  double
    r = this->Udouble(0, 1),
    Praise = pow(1 + 0.5*sqrt(M_PI)/this->alpha, -1.0);

  if (Praise > r) {
    this->raiseOrder();
  } else {
    this->lowerOrder();
  }




}