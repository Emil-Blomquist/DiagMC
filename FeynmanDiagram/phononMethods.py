import sys

sys.path.append('../')
from Propagator import D

def addInternalPhonon(self, v1, v2, momentum, sinTheta):
  d = D(momentum, sinTheta)

  if v1.position > v2.position:
    print('ERROR: Phonon cannot run backwards')
    return False

  d.setStart(v1)
  d.setEnd(v2)

  # append to list
  self.Ds.append(d)

  # change momentum for electrons under phonon for conservation of momentum
  g = v1.G[1]
  while g.start != v2:
    g.addMomentum(-momentum)
    g = g.end.G[1]

  # return phonon
  return d

def removeInternalPhonon(self, dIndex):
  d = self.Ds[dIndex]

  # change momentum for electrons under phonon for conservation of momentum
  g = d.start.G[1]
  while g.start != d.end:
    g.addMomentum(d.momentum)
    g = g.end.G[1]

  # will unlink from vertices
  d.remove()

  # remove from Ds list                                                            <-- takes time => hash table
  del self.Ds[dIndex]

def setInternalPhononMomentum(self, d, momentum, sinTheta):
  # sinTheta needed for Jacobian
  d.setSinTheta(sinTheta)

  dP = momentum - d.momentum

  # add momentum difference
  d.addMomentum(dP)

  # change momentum for electrons under phonon for conservation of momentum
  g = d.start.G[1]
  while g.start != d.end:
    g.addMomentum(-dP)
    g = g.end.G[1]

  # return changed phonon
  return d

def swapPhononEnds(self, v1, v2):
  # order temporally
  if v1.position > v2.position:
    (v1, v2) = (v2, v1)

  # need to be neighbours
  if v1.G[1].end != v2:
    print('ERROR: vertices not neighbours')
    return False
    
  d1 = v1.D[0] or v1.D[1]
  d2 = v2.D[0] or v2.D[1]

  # change in momentum
  dp1 = -d1.momentum if v1.D[0] else d1.momentum
  dp2 = -d2.momentum if v2.D[1] else d2.momentum
  dp = dp1 + dp2

  # conserve momentum
  v1.G[1].addMomentum(dp)

  # relink
  if d1 == d2:
    # just reverse the momentum
    d1.setMomentum(-d1.momentum)
  else:
    if d1.start == v1:
      d1.setStart(v2)
    else:
      d1.setEnd(v2)

    if d2.start == v2:
      d2.setStart(v1)
    else:
      d2.setEnd(v1)