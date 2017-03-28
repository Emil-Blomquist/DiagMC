#include "../DiagrammaticMonteCarlo.h"

void DiagrammaticMonteCarloV2::changeDiagramLengthComplex (double param) {
  // pick a electron propagator on random
  shared_ptr<Electron> g = this->FD.Gs[this->Uint(0, this->FD.Gs.size() - 1)];

  

  
}