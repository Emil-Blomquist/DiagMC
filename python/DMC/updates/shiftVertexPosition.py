import numpy as np

def shiftVertexPosition(self):
  # pick a vertex on random
  g = np.random.choice(self.FD.Gs[0:-1])
  v = g.end

  # fetch available time interval
  t1 = g.start.position
  t2 = v.G[1].end.position

  c = -1 if v.D[0] else 1
  dE = 0.5*np.linalg.norm(v.G[0].momentum)**2 - 0.5*np.linalg.norm(v.G[1].momentum)**2 - c
  dtdE = (t2 - t1)*dE

  # sample t
  r = np.random.rand()
  if -dtdE > 100:
    # to avoid overflow due to exponential
    t = t2 - np.log(r)/dE
  else:
    t = t1 - np.log(1 - r*(1 - np.exp(-dtdE)))/dE

  if self.debug:
    tOld = v.position
    diagOld = self.FD()

  self.FD.setVertexPosition(v, t)

  if self.debug:
    diag = self.FD()

    R = np.exp(np.log(diag) - np.log(diagOld) + dE*(t - tOld))

    print('Shift Vertex Position', R)

  return 1