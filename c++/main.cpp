#include "FeynmanDiagram/FeynmanDiagram.h"

int main() {
  FeynmanDiagram FD(Vector3d(0,0,0), 1, 2, 3);

  FD.insertVertex(FD.Vs.begin(), 0.1);
  FD.insertVertex(next(FD.Vs.begin(), 1), 0.1);

  FD.addInternalPhonon(next(FD.Vs.begin(), 1), prev(FD.Vs.end(), 2), Vector3d(1,2,1), 1, 2);
  FD.addInternalPhonon(FD.Vs.begin(), prev(FD.Vs.end(), 1), Vector3d(3,4,5), 1, 2);

  // FD.print();

  FD.plot();
}