#include <iostream>
#include <eigen3/Eigen/Dense>

#include "FeynmanDiagram/FeynmanDiagram.h"
#include "Vertex/Vertex.h"

int main() {
  FeynmanDiagram FD(Vector3d(1,2,3), 1, 2, 3);


  FD.start->print();
  FD.end->print();
}