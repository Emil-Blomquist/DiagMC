import numpy as np

def calculateP0(d):
  dt = d.end.position - d.start.position

  P0 = np.array([0, 0, 0])

  # loop through electrons between phonon propagator ends
  g = d.start.G[1]
  while g.start != d.end:
    P0 = P0 + g.momentum*(g.end.position - g.start.position)
    g = g.end.G[1]

  # add own momentum as well
  P0 = P0/dt + d.momentum

  return P0