import numpy as np

class Propagator(object):
  
  def __init__(self, which, momentum):
    self.type = which
    self.momentum = momentum
    self.start = None
    self.end = None

    # used to include Jacobian
    self.sinTheta = None

    # saved variables
    self.momentumSaved = None
    self.startSaved = None
    self.endSaved = None

  def setStart(self, vertex):
    # unlink from previous vertex
    if self.start and getattr(self.start, self.type)[1] == self:
      getattr(self.start, self.type)[1] = None

    self.start = vertex
    vertex.setPropagatorOut(self)

  def setEnd(self, vertex):
    # unlink from previous vertex
    if self.end and getattr(self.end, self.type)[0] == self:
      getattr(self.end, self.type)[0] = None

    self.end = vertex
    vertex.setPropagatorIn(self)

  def addMomentum(self, p):
    self.momentum = self.momentum + p

  def setMomentum(self, p):
    self.momentum = p

  def remove(self):
    # unlink from vertices
    getattr(self.start, self.type)[1] = None
    getattr(self.end, self.type)[0] = None

  def save(self):
    self.momentumSaved = self.momentum[:]
    self.startSaved = self.start
    self.endSaved = self.end
    self.sinThetaSaved = self.sinTheta

  def revert(self):
    self.momentum = self.momentumSaved
    self.start = self.startSaved
    self.end = self.endSaved
    self.sinTheta = self.sinThetaSaved

  def __call__(self):
    k2 = np.linalg.norm(self.momentum)**2
    t = self.end.position - self.start.position

    if t < 0:
      print('ERROR: propagator running backwards')
      return False

    if self.type == 'G':
      return np.exp(-0.5*k2*t)
    else:
      # coupling constant
      a = 2

      print('call', self.momentum, self.sinTheta)

      if True:
        # internal phonon

        # Jacobian and vertex contribution
        # omit minus sign here since it will cancel with another one when summing diagrams together
        contr = 8**0.5 * a * np.pi * self.sinTheta

        return contr * np.exp(-t)









class G(Propagator):
  def __init__(self, momentum):
    Propagator.__init__(self, 'G', momentum)

class D(Propagator):
  def __init__(self, momentum):
    Propagator.__init__(self, 'D', momentum)

  def getP0(self):
    dt = self.end.position - self.start.position

    P0 = np.array([0, 0, 0])

    # loop through electrons under phonon arc
    g = self.start.G[1]
    while g.start != self.end:
      P0 = P0 + g.momentum*(g.end.position - g.start.position)
      g = g.end.G[1]

    # add own momentum as well
    P0 = P0/dt + self.momentum

    return P0

  def recalculateSinTheta(self):
    P0 = self.getP0()
    p0 = np.linalg.norm(P0)
    q = np.linalg.norm(self.momentum)

    # if P0 = 0 we use that theta is the angle against the z-axis
    # we need 10^-10 to handle rounding errors
    if p0 > 10**-10:
      print('p0>0', p0)
      print('WAT?', np.dot(self.momentum, P0)/(q*p0))
      cosTheta = np.dot(self.momentum, P0)/(q*p0)
    else:
      print('p0=0', p0)
      cosTheta = self.momentum[2]/q

    # rounding error
    if cosTheta**2 > 1:
      sinTheta = 0
    else:
      sinTheta = (1 - cosTheta**2)**0.5

    self.sinTheta = sinTheta

    # need to return
    return sinTheta

    


