import numpy as np
from math import factorial

def lowerOrder(self):
  # save old diagram value
  diagOld = self.FD()

  # choose phonon on random
  dIndex = np.random.randint(0, len(self.FD.Ds))
  d = self.FD.Ds[dIndex]

  # inverse probabilities for internal parameters
  wInvt1 = d.start.G[1].end.position - d.start.G[0].start.position

  dt = self.FD.time - d.start.position
  t = d.end.position - d.start.position
  wInvt = (1 - np.exp(-dt))*np.exp(t)

  std = t**-0.5;
  q = np.linalg.norm(d.momentum)
  # wInvQ = 2*np.pi**2 * (np.pi/2*std**2)**0.5*np.exp(0.5*(q/std)**2)

  # temporarily save vertices
  v1 = d.start
  v2 = d.end

  # remove phonon
  self.FD.removeInternalPhonon(dIndex)

  # remove vertices
  self.FD.removeVertex(v1)
  self.FD.removeVertex(v2)


  # save diagarm order
  n = len(self.FD.Ds)

  # ration
  r = (n + 1)/factorial(2*n + 1)

  # get current diagram value
  diag = self.FD()

  # R = r*diag/diagOld / (wInvt1 * wInvt * wInvQ)
  R = r

  return R