import numpy as np

class Propagator(object):
  
  def __init__(self, which, momentum, sinTheta = None):
    self.type = which
    self.momentum = momentum
    self.start = None
    self.end = None

    # used for jacobian, for phonons only
    self.sinTheta = sinTheta

    # saved variables
    self.momentumSaved = None
    self.sinThetaSaved = None
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
    self.sinThetaSaved = self.sinTheta
    self.startSaved = self.start
    self.endSaved = self.end

  def revert(self):
    self.momentum = self.momentumSaved
    self.sinTheta = self.sinThetaSaved
    self.start = self.startSaved
    self.end = self.endSaved

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
  def __init__(self, momentum, sinTheta):
    Propagator.__init__(self, 'D', momentum, sinTheta)

  def setSinTheta(self, sinTheta):
    self.sinTheta = sinTheta

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

  def sinTheta(self):
    ##
    ## if P0 = 0 we use that theta is the angle against the z-axis
    ##


