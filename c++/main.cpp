#include <iostream>

#include "Propagator/G.h"
#include "Vertex/Vertex.h"
 
int main()
{
  G g1(Vector3d(1,2,3), 1.1, 1.2), g2(Vector3d(3,2,1), 2.1, 2.2);
  Vertex v1(1), v2(2);

  g1.setStart(&v1);
  g1.setEnd(&v2);

  v2.setIngoingG(&g1);
  v2.setOutgoingG(&g2);

  v1.print();
  // v2.print();
  g1.print();
  g2.print();

}