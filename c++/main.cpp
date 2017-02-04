#include "FeynmanDiagram/FeynmanDiagram.h"

int main() {
  FeynmanDiagram FD(Vector3d(0,0,0), 1, 2, 3);

  FD.insertVertex(FD.Vs.begin(), 0.1);
  FD.addInternalPhonon(FD.Vs.begin(), --FD.Vs.end(), Vector3d(1,2,3), 1, 2);

  // FD.print();

  FD.plot();
}