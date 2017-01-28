import sys
import numpy as np

sys.path.append('/dependencies')
from .dependencies.calculateQ import calculateQ
from .dependencies.calculateP0 import calculateP0

def changeInternalPhononMomentumDirection(self):
  # choose phonon propagator on random
  d = np.random.choice(self.FD.Ds)

  q = np.linalg.norm(d.momentum)
  P0 = calculateP0(d)
  p0 = np.linalg.norm(P0)
  
  r = np.random.rand()
  phi = np.random.uniform(0, 2*np.pi)

  ##
  ## sample new theta
  ##
  a = q*p0*(d.end.position - d.start.position)
  if a < 10**-10:
    # small a approximation of the one below
    cosTheta = 1 - 2*r
  else:
    cosTheta = 1 + np.log(1 - r*(1 - np.exp(-2*a)))/a
  theta = np.arccos(cosTheta)

  # generate new Q for phonon momentum
  Q = calculateQ(P0, q, theta, phi)

  if self.debug:

    # if P0 = 0 we use that theta is the angle against the z-axis
    # we need 10^-10 to handle rounding errors
    if p0 > 10**-10:
      cosThetaOld = np.dot(d.momentum, P0)/(q*p0)
    else:
      cosThetaOld = d.momentum[2]/q

    if cosThetaOld <= -1:
      thetaOld = np.pi
    else:
      thetaOld = np.arccos(cosThetaOld)

    ##
    ## "d.theta" is outdated since the diagram might have changed in the sense that
    ## it is no longer valid for comparing Q with a new Q'. Even the same Q might be
    ## more/less preffered. "d.theta" is only good for evaluating the diagram.
    ##

    ##
    ## Q no longer corresponds to theta since P0 might have changed:                        TRY TO CHANGE ALL THETAS AFTER EACH DIAGRAM UPDATE
    ##                                                                                      -> might give better data
    ##                                                                                      -> we cannot, then we need to take that into account to get R = 1
    ##
    self.FD.setInternalPhononMomentumAngle(d, thetaOld, d.phi)

    # calculate the value correcponding to the old diagram
    diagOld = self.FD()

  # update diagram
  self.FD.setInternalPhononMomentum(d, Q)
  self.FD.setInternalPhononMomentumAngle(d, theta, phi)

  if self.debug:
    diag = self.FD()

    if abs(diagOld) > 0 and np.sin(theta) > 0:
      R = diag/diagOld
      R *= np.sin(thetaOld)/np.sin(theta)
      R *= np.exp(-a*(np.cos(theta) - np.cos(thetaOld))) 
      print('Momentum Direction', R)
    else:
      R = 1
      print('--------------------------------------------------------------------')
      print('Momentum Direction: overflow', R)
      print('--------------------------------------------------------------------')

  return 1