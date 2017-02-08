#include "FeynmanDiagram/FeynmanDiagram.h"
#include "Display/Display.h"

int main() {
  FeynmanDiagram FD(Vector3d(0,0,0), 1, 2, 3);

  auto v1 = FD.insertVertex(0, 0.1);
  auto v2 = FD.insertVertex(1, 0.1);
  auto v3 = FD.insertVertex(2, 0.1);
  auto v4 = FD.insertVertex(3, 0.1);

  auto d1 = FD.addInternalPhonon(v1, v2, Vector3d(1,2,1), 1, 2);
  auto d2 = FD.addInternalPhonon(v3, v4, Vector3d(3,4,5), 1, 2);


  // FD.setVertexPosition(v2, 0.15);

  // FD.removeInternalPhonon(d1);
  // FD.removeVertex(v1);
  // FD.removeVertex(v3);

  // FD.setInternalPhononMomentum(d1, Vector3d(10,10,10));
  // FD.setInternalPhononMomentumDirection(d1, 3, 2);

  // FD.swapPhonons(v1, v2);
  FD.swapPhonons(v2, v3);
  FD.swapPhonons(v1, v2);
  FD.swapPhonons(v3, v4);
  FD.swapPhonons(v2, v3);


  Display display(&FD);
  display.render();
}