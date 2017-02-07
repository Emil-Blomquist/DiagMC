#include "FeynmanDiagram/FeynmanDiagram.h"
#include "plot/plot.h"

int main() {
  FeynmanDiagram FD(Vector3d(0,0,0), 1, 2, 3);

  auto v1 = FD.insertVertex(0, 0.1);
  auto v2 = FD.insertVertex(1, 0.1);
  auto v3 = FD.insertVertex(2, 0.1);
  auto v4 = FD.insertVertex(3, 0.1);

  auto d1 = FD.addInternalPhonon(v1, v3, Vector3d(1,2,1), 1, 2);
  auto d2 = FD.addInternalPhonon(v2, v4, Vector3d(3,4,5), 1, 2);

  plot plt(&FD);
  plt.render();

}