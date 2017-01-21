import numpy as np

def raiseOrder(self):
  # old diagram value
  diagOld = self.FD()

  # select electron lines to split
  if len(self.FD.Ds) == 0:
    i1, i2 = 0, 1
  else:
    i1, i2 = 0, 2

  # split first electron line
  g1 = self.FD.Gs[i1]
  t1 = np.random.uniform(g1.start.position, g1.end.position)
  wInvt1 = g1.end.position - g1.start.position

  dt1 = t1 - g1.start.position
  v1 = self.FD.insertVertex(i1, dt1)

  # split second electron line
  g2 = self.FD.Gs[i2]
  t2 = np.random.uniform(g2.start.position, g2.end.position)
  wInvt2 = g2.end.position - g2.start.position

  dt2 = t2 - g2.start.position
  v2 = self.FD.insertVertex(i2, dt2)

  # sample phonon momentum
  std = (t2 - t1)**-0.5;
  q = np.abs(np.random.normal(0, std))
  theta = np.random.uniform(0, np.pi)
  phi = np.random.uniform(0, 2*np.pi)
  wInvQ = 2*np.pi**2 * (np.pi*std**2 / 2)**0.5 * np.exp(0.5*(q/std)**2)

  Q = np.array([
    q*np.sin(theta)*np.cos(phi),
    q*np.sin(theta)*np.sin(phi),
    q*np.cos(theta)
  ])

  # add phonon
  self.FD.addInternalPhonon(v1, v2, Q, theta, phi)

  # current diagram value
  diag = self.FD()

  # acceptance ratio
  R = diag/diagOld * (wInvt1 * wInvt2 * wInvQ)

  return R