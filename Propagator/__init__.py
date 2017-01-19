import numpy as np

class Propagator(object):
  
  def __init__(self, which, momentum):
    self.type = which
    self.momentum = momentum
    self.start = None
    self.end = None

    # used to include Jacobian
    self.theta = None
    self.phi = None

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

  def unlink(self):
    # unlink from vertices
    getattr(self.start, self.type)[1] = None
    getattr(self.end, self.type)[0] = None

  def save(self):
    self.momentumSaved = self.momentum[:]
    self.startSaved = self.start
    self.endSaved = self.end
    self.thetaSaved = self.theta
    self.phiSaved = self.phi

  def revert(self):
    self.momentum = self.momentumSaved
    self.start = self.startSaved
    self.end = self.endSaved
    self.theta = self.thetaSaved
    self.phi = self.phiSaved

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
        contr = a/(8**0.5*np.pi**2) * np.sin(self.theta)

        return contr * np.exp(-t)









class G(Propagator):
  def __init__(self, momentum):
    Propagator.__init__(self, 'G', momentum)

class D(Propagator):
  def __init__(self, momentum):
    Propagator.__init__(self, 'D', momentum)

  def setTheta(self, theta):
    self.theta = theta

  def setPhi(self, phi):
    self.phi = phi

    


