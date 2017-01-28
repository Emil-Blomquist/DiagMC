from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt

from scipy.special import erfinv, erf

import warnings
warnings.filterwarnings('error')

sys.path.append('../')

from FeynmanDiagram import FeynmanDiagram

# print(sys.float_info.max, np.log(sys.float_info.max))


##
##
## If we use the proper probability distribution when sampling
## internal variables we wont need to save the invProb in the
## feynmanDiagram object since all of that will be taken care of here
##
##

##
## swap phonon ends
##
def swapPhononEnds(fd):
  # pick an internal electron propagator on random
  g = np.random.choice(fd.Gs[1:-1])

  v1 = g.start
  v2 = g.end

  c1 = 1 if v1.D[1] else -1
  c2 = 1 if v2.D[1] else -1

  d1 = v1.D[0] or v1.D[1]
  d2 = v2.D[0] or v2.D[1]

  t = v2.position - v1.position
  Eafter = 0.5*np.linalg.norm((g.momentum + c1*d1.momentum - c2*d2.momentum))**2
  Ebefore = 0.5*np.linalg.norm(g.momentum)**2

  # if both vertices belong to the same phonon
  # same amount of phonon is still present above after the swap in this case
  if d1 == d2:
    dw = 0
  else:
    dw = c2 - c1

  try:
    R = np.exp(-t*(Eafter - Ebefore + dw))
  except Warning:
    R = 1
    print('ERROR: swapPhononEnds', -t*(Eafter - Ebefore + dw))

  fd.swapPhononEnds(v1, v2)

  return R


##
## change phonon momentum
##
def changePhononMomentumMagnitudeNaive(fd):
  # pick an internal phonon propagator on random
  d = np.random.choice(fd.Ds)

  standDev = (d.end.position - d.start.position)**-0.5;
  q = np.abs(np.random.normal(0, standDev));

  Q = np.array([
    q*np.sin(d.theta)*np.cos(d.phi),
    q*np.sin(d.theta)*np.sin(d.phi),
    q*np.cos(d.theta)
  ]);

  diagOld = fd()
  qOld = np.linalg.norm(d.momentum)


  fd.setInternalPhononMomentum(d, Q)

  diag = fd()

  if diagOld == 0 and diag == 0:
    R = 0.5
  elif diagOld == 0:
    R = 1
  else:
    R = diag/diagOld * np.exp(0.5*(q/standDev)**2 - 0.5*(qOld/standDev)**2)

  return R

def changePhononMomentumDirectionNaive(fd):
  # pick an internal phonon propagator on random
  d = np.random.choice(fd.Ds)

  q = np.linalg.norm(d.momentum)
  theta = np.random.uniform(0, np.pi);
  phi = np.random.uniform(0, 2*np.pi);

  Q = np.array([
    q*np.sin(theta)*np.cos(phi),
    q*np.sin(theta)*np.sin(phi),
    q*np.cos(theta)
  ]);



  diagOld = fd()

  fd.setInternalPhononMomentumAngle(d, theta, phi)
  fd.setInternalPhononMomentum(d, Q)

  diag = fd()


  if diagOld == 0 and diag == 0:
    R = 0.5
  elif diagOld == 0:
    R = 1
  else:
    R = diag/diagOld

  return R




def shiftVertexPosition(fd):
  # pick a vertex on random
  g = np.random.choice(fd.Gs[0:-1])
  v = g.end

  # fetch available time interval
  t1 = g.start.position
  t2 = v.G[1].end.position

  c = -1 if v.D[0] else 1
  dE = 0.5*np.linalg.norm(v.G[0].momentum)**2 - 0.5*np.linalg.norm(v.G[1].momentum)**2 - c

  # approximation (good one) to solve "RuntimeWarning: overflow encountered in exp"
  r = np.random.rand()
  if -(t2 - t1)*dE > 100:
    t = t2 - np.log(r)/dE
  else:
    t = t1 - np.log(1 - r*(1 - np.exp(-(t2 - t1)*dE)))/dE

  fd.setVertexPosition(v, t)

  return 1

def calculateP0(d):
  dt = d.end.position - d.start.position

  P0 = np.array([0, 0, 0])

  # loop through electrons under phonon arc
  g = d.start.G[1]
  while g.start != d.end:
    P0 = P0 + g.momentum*(g.end.position - g.start.position)
    g = g.end.G[1]

  # add own momentum as well
  P0 = P0/dt + d.momentum

  return P0

num = 0
def QofP0(P0, q, theta, phi):
  global num

  p0 = np.linalg.norm(P0)

  # 10^-10 to handle rounding errors
  # if p0 > 10**-14:
  #   Ep = P0/p0
  # else:
  #   num += 1
  #   # since P0 has no length we choose this as the direction
  #   Ep = np.array([0, 0, 1])

  try:
    Ep = P0/p0
  except Warning:
    num += 1
    Ep = np.array([0, 0, 1])


  # generate a temporary vector used to obtain two other vectors in order to span the rest of R^3
  tempVector = Ep + np.array([1, 0, 0])
  Eo1 = np.cross(Ep, tempVector)
  if np.linalg.norm(Eo1) < 10**-14:
    tempVector = Ep + np.array([0, 1, 0])
    Eo1 = np.cross(Ep, tempVector)

  Eo1 = Eo1/np.linalg.norm(Eo1)
  Eo2 = np.cross(Ep, Eo1)

  Qp = Ep * q*np.cos(theta)
  Qo = (Eo1*np.cos(phi) + Eo2*np.sin(phi)) * q*np.sin(theta)

  return Qp + Qo

def changePhononMomentumDirection(fd):
  global debug

  # choose phonon propagator on random
  d = np.random.choice(fd.Ds)

  q = np.linalg.norm(d.momentum)
  P0 = calculateP0(d)
  p0 = np.linalg.norm(P0)
  
  r = np.random.rand()
  phi = np.random.uniform(0, 2*np.pi)

  # 10^-10 to handle rounding errors
  
  a = q*p0*(d.end.position - d.start.position)
  if a < 10**-10:
    cosTheta = 1 - 2*r
  else:
    cosTheta = 1 + np.log(1 - r*(1 - np.exp(-2*a)))/a
  theta = np.arccos(cosTheta)

  Q = QofP0(P0, q, theta, phi)

  if debug:
    a = q*p0*(d.end.position - d.start.position)

    # if P0 = 0 we use that theta is the angle against the z-axis
    # we need 10^-10 to handle rounding errors
    if p0 > 10**-14:
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
    fd.setInternalPhononMomentumAngle(d, thetaOld, d.phi)

    diagOld = fd()

  # update diagram
  fd.setInternalPhononMomentum(d, Q)
  fd.setInternalPhononMomentumAngle(d, theta, phi)

  if debug:
    diag = fd()
    a = q*p0*(d.end.position - d.start.position)

    if abs(diagOld) > 0 and np.sin(theta) > 0:
      R = diag/diagOld
      R *= np.sin(thetaOld)/np.sin(theta)
      R *= np.exp(-a*(np.cos(theta) - np.cos(thetaOld))) 
    else:
      R = 1

    # print('Momentum Direction', R)

  return 1

def changePhononMomentumMagnitud(fd):
  global debug

  # choose phonon propagator on random
  d = np.random.choice(fd.Ds)

  P0 = calculateP0(d)
  p0 = np.linalg.norm(P0)

  sqrta = (0.5*(d.end.position - d.start.position))**0.5
  b = p0*np.cos(d.theta)
  r = np.random.rand()
  q = b + erfinv(r + (r - 1)*erf(sqrta*b))/sqrta
  Q = QofP0(P0, q, d.theta, d.phi)

  if debug:
    qOld = np.linalg.norm(d.momentum)
    QOld = QofP0(P0, qOld, d.theta, d.phi)
    fd.setInternalPhononMomentum(d, QOld)
    diagOld = fd()

  # update diagram
  fd.setInternalPhononMomentum(d, Q)
  diag = fd()

  if debug:

    try:
      R = diag/diagOld
      R *= np.exp(-sqrta**2 * ((qOld - b)**2 - (q - b)**2))
    except Warning:
      print(diag, diagOld, np.exp(-sqrta**2 * ((qOld - b)**2 - (q - b)**2)))


    # print('Momentum Magnitude', R)

  return 1







##
##
##







def DMC(t, P0, N, sophisticated = True):

  # create second order diagram
  feynmanDiagram = FeynmanDiagram(t, P0)

  dt = t/5
  v1 = feynmanDiagram.insertVertex(0, dt)
  v2 = feynmanDiagram.insertVertex(1, dt)
  v3 = feynmanDiagram.insertVertex(2, dt)
  v4 = feynmanDiagram.insertVertex(3, dt)
  feynmanDiagram.addInternalPhonon(v1, v3, P0, 1, 1) 
  feynmanDiagram.addInternalPhonon(v2, v4, P0, 1, 1)

  # to reach some sense of randomness
  # for i in range(0, 10):
  #   if sophisticated:
  #     changePhononMomentumDirection(feynmanDiagram)
  #     changePhononMomentumMagnitud(feynmanDiagram)
  #   else:
  #     changePhononMomentumDirectionNaive(feynmanDiagram)
  #     changePhononMomentumMagnitudeNaive(feynmanDiagram)

  #   shiftVertexPosition(feynmanDiagram)
  #   swapPhononEnds(feynmanDiagram)

  if sophisticated:
    updates = [shiftVertexPosition, swapPhononEnds, changePhononMomentumMagnitud, changePhononMomentumDirection]
  else:
    updates = [shiftVertexPosition, swapPhononEnds, changePhononMomentumMagnitudeNaive, changePhononMomentumDirectionNaive]

  # updates = [swapPhononEnds]

  Rnaive = []
  Rswap = []

  bins = {}
  for i in range(0, N):
    feynmanDiagram.save()

    update = np.random.choice(updates)

    r = update(feynmanDiagram)

    if update == changePhononMomentumMagnitudeNaive or update == changePhononMomentumDirectionNaive:
      Rnaive.append(r)
    elif update == swapPhononEnds:
      Rswap.append(r)

    if np.random.rand() > min(1, r):
      # reject step
      feynmanDiagram.revert()

    struct = feynmanDiagram.structure()
    if struct in bins:
      bins[struct] += 1
    else:
      bins[struct] = 1

  for key in bins:
    bins[key] = bins[key]/N

  return bins, Rnaive, Rswap









##
##
##








debug = False


P0 = np.array([0, 0, 6])

n = 20
N = 1000

a1122 = []
a1212 = []
a1221 = []
b1122 = []
b1212 = []
b1221 = []

plt.ion()

for i in range(0, n):
  bins = DMC(1, P0, N, True)

  a1122.append(bins[1122])
  a1212.append(bins[1212])
  a1221.append(bins[1221])

  plt.plot(a1122, 'b-')
  plt.plot(a1212, 'b-')
  plt.plot(a1221, 'b-')
  plt.pause(0.05)

print([a1122, a1212, a1221])


for i in range(0, n):
  bins = DMC(1, P0, N, False)

  b1122.append(bins[1122]) if 1122 in bins else b1122.append(0)
  b1212.append(bins[1212]) if 1212 in bins else b1212.append(0)
  b1221.append(bins[1221]) if 1221 in bins else b1221.append(0)

  plt.plot(b1122, 'r-')
  plt.plot(b1212, 'r-')
  plt.plot(b1221, 'r-')
  plt.pause(0.05)

print([b1122, b1212, b1221])
print(num)

plt.ioff()

plt.plot([0, n], [np.mean(a1122)] * 2, 'b--')
plt.plot([0, n], [np.mean(a1212)] * 2, 'b--')
plt.plot([0, n], [np.mean(a1221)] * 2, 'b--')
plt.plot([0, n], [np.mean(b1122)] * 2, 'r--')
plt.plot([0, n], [np.mean(b1212)] * 2, 'r--')
plt.plot([0, n], [np.mean(b1221)] * 2, 'r--')


plt.plot(a1122, 'b-', label='a1122')
plt.plot(a1212, 'b-', label='a1212')
plt.plot(a1221, 'b-', label='a1221')
plt.plot(b1122, 'r-', label='b1122')
plt.plot(b1212, 'r-', label='b1212')
plt.plot(b1221, 'r-', label='b1221')

print(np.mean(a1122), np.mean(b1122))
print(np.std(a1122), np.std(b1122))
print(np.std(a1212), np.std(b1212))
print(np.std(a1221), np.std(b1221))



plt.savefig('plot.pdf')

plt.show()