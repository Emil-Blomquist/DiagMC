# import sys
import numpy as np
# from math import factorial

# sys.path.append('/dependencies')
# from .dependencies.numDiagramsWithPhonon import numDiagramsWithPhonon
# from .dependencies.numDiagramsWithoutPhonon import numDiagramsWithoutPhonon

def lowerOrder(self):
  # save old diagram value
  diagOld = self.FD()
  # save number of diagrams we might go to by removing a phonon
  # nDown = numDiagramsWithoutPhonon(self.FD)

  # choose phonon on random
  # dIndex = np.random.randint(0, len(self.FD.Ds))
  # d = self.FD.Ds[dIndex]

  n = len(self.FD.Ds)

  # requirement for lowering: first and third vertex belong to the same phonon
  if n > 1 and self.FD.Gs[0].end.D[1] != self.FD.Gs[2].end.D[0]:
    return 0


  d = self.FD.Gs[0].end.D[1]




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
  self.FD.removeInternalPhonon(d)

  # remove vertices
  self.FD.removeVertex(v1)
  self.FD.removeVertex(v2)

  # diagarm order and how many diagrams we might go to
  # n = len(self.FD.Ds)
  # nUpp = numDiagramsWithPhonon(self.FD)

  # # ration
  # if n == 0:
  #   r = 1
  # else:
  #   r = (n + 1)/(2*n) * nDown/nUpp

  # get current diagram value
  diag = self.FD()

  r = 1

  # R = r*diag/diagOld / (wInvt1 * wInvt * wInvQ)
  R = r

  return R