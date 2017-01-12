import sys
import numpy as np

sys.path.append('../')
from Vertex import Vertex
from Propagator import G

def insertVertex(self, g2splitIndex, dt):
  g1 = self.Gs[g2splitIndex]
  g2 = G(np.copy(g1.momentum))

  if dt < 0 or dt > g1.end.position - g1.start.position:
    print('ERROR: Propagator not present at splitting time')
    return False

  # relink
  g2.setEnd(g1.end)
  g1.setEnd(Vertex(g1.start.position + dt))
  g2.setStart(g1.end)

  # insert new propagator into chronologically ordered list
  self.Gs.insert(g2splitIndex + 1, g2)

  # returns the inserted vertex
  return g1.end

def removeVertex(self, vertexIndex):
  v = self.Gs[vertexIndex].end

  if v.D[0] or v.D[1]:
    print('ERROR: Phonon line still attatched to vertex')
    return False

  if np.linalg.norm(v.G[1].momentum - v.G[0].momentum) != 0:
    print('ERROR: Electron lines do not conserve momentum over vertex')
    return False

  # relink
  v.G[0].setEnd(v.G[1].end)

  # remove second G from Gs
  del self.Gs[vertexIndex + 1]

def setVertexPosition(self, v, position):
  if not v.G[0].start.position < position < v.G[1].end.position:
    print('ERROR: New vertex time invalid')
    return False

  v.setPosition(position)

  # return changed vertex
  return v