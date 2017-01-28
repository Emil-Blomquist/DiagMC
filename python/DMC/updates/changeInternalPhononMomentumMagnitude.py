import sys
import numpy as np
from scipy.special import erfinv, erf

sys.path.append('/dependencies')
from .dependencies.calculateQ import calculateQ
from .dependencies.calculateP0 import calculateP0

def changeInternalPhononMomentumMagnitude(self):
  # choose phonon propagator on random
  d = np.random.choice(self.FD.Ds)

  P0 = calculateP0(d)
  p0 = np.linalg.norm(P0)

  sqrta = (0.5*(d.end.position - d.start.position))**0.5
  b = p0*np.cos(d.theta)
  r = np.random.rand()
  q = b + erfinv(r + (r - 1)*erf(sqrta*b))/sqrta

  Q = calculateQ(P0, q, d.theta, d.phi)

  if self.debug:
    qOld = np.linalg.norm(d.momentum)
    QOld = calculateQ(P0, qOld, d.theta, d.phi)
    self.FD.setInternalPhononMomentum(d, QOld)
    diagOld = self.FD()

  # update diagram
  self.FD.setInternalPhononMomentum(d, Q)

  if self.debug:
    diag = self.FD()

    R = diag/diagOld
    R *= np.exp(-sqrta**2 * ((qOld - b)**2 - (q - b)**2))

    print('Momentum Magnitude', R)

  return 1