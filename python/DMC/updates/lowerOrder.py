import numpy as np

def lowerOrder(self):
  # first encountered phonon line
  d = self.FD.Gs[0].end.D[1]

  # requirement to lower: first and third vertex belong to the same phonon
  if len(self.FD.Ds) > 1 and d != self.FD.Gs[2].end.D[0]:
    return 0

  # old diagram value
  diagOld = self.FD()

  # vertices
  v1 = d.start
  v2 = d.end

  # inverse probabilities for internal parameters
  if len(self.FD.Ds) == 1:
    wInvt1 = self.FD.time
    wInvt2 = d.end.G[1].end.position - d.end.G[0].start.position
  else:
    wInvt1 = d.start.G[1].end.position - d.start.G[0].start.position
    wInvt2 = d.end.G[1].end.position - d.end.G[0].start.position

        # dt = self.FD.time - d.start.position
        # t = d.end.position - d.start.position
        # wInvt = (1 - np.exp(-dt))*np.exp(t)

  std = (v2.position - v1.position)**-0.5;
  q = np.linalg.norm(d.momentum)
  wInvQ = 2*np.pi**2 * (0.5*np.pi*std**2)**0.5 * np.exp(0.5*(q/std)**2)

  # remove phonon
  self.FD.removeInternalPhonon(d)

  # remove vertices
  self.FD.removeVertex(v1)
  self.FD.removeVertex(v2)

  # get current diagram value
  diag = self.FD()

  # acceptance ratio
  # if diagOld == 0 or wInvQ == 0:
  #   R = 1
  # else:
  R = diag/diagOld / (wInvt1 * wInvt2 * wInvQ)

  return R