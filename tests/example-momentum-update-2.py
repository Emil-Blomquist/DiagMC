from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

from scipy.special import erfinv, erf

import warnings
warnings.filterwarnings('error')

from FeynmanDiagram import FeynmanDiagram

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
# def swapPhononEnds(fd):
#   # pick an internal electron propagator on random
#   g = np.random.choice(fd.Gs[1:-1])

#   v1 = g.start
#   v2 = g.end

#   c1 = 1 if v1.D[1] else -1
#   c2 = 1 if v2.D[1] else -1

#   d1 = v1.D[0] or v1.D[1]
#   d2 = v2.D[0] or v2.D[1]

#   t = v2.position - v1.position
#   Eafter = 0.5*np.linalg.norm((g.momentum + c1*d1.momentum - c2*d2.momentum))**2
#   Ebefore = 0.5*np.linalg.norm(g.momentum)**2

#   # if both vertices belong to the same phonon
#   # same amount of phonon is still present above after the swap in this case
#   if d1 == d2:
#     dw = 0
#   else:
#     dw = c2 - c1

#   try:
#     R = np.exp(-t*(Eafter - Ebefore + dw))
#   except Warning:
#     R = 1
#     print('ERROR: swapPhononEnds', -t*(Eafter - Ebefore + dw))

#   fd.swapPhononEnds(v1, v2)

#   return R


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

  # R = np.exp(0.5*(q/standDev)**2)*diag / (np.exp(0.5*(qOld/standDev)**2)*diagOld)

  ##
  ## build try catch:
  ## RuntimeWarning: divide by zero encountered in double_scalars
  ## RuntimeWarning: invalid value encountered in double_scalars
  ##

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

    print('Momentum Direction', R)

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


    print('Momentum Magnitude', R)

  return 1







##
##
##







def DMC(fd, N, sophisticated = True):
  if sophisticated:
    updates = [changePhononMomentumMagnitud, changePhononMomentumDirection]
  else:
    updates = [changePhononMomentumMagnitudeNaive, changePhononMomentumDirectionNaive]

  D = []
  M = []
  T = []
  P = []

  feynmanDiagram.save()
  for i in range(0, N):


    update = np.random.choice(updates)

    
    r = update(fd)

    if np.random.rand() > min(1, r):
      # reject step
      fd.revert()
    else:
      fd.save()

    D.append(fd())
    M.append(np.linalg.norm(fd.Ds[0].momentum))
    T.append(fd.Ds[0].theta)
    P.append(fd.Ds[0].phi)

  return D, M, T, P



##
## divide by W and then append to count
##







debug = False

t = 1
P = np.array([1, 0, 1])
dt = t/3

nNaive = 100000
nSoph = 10000




# create the figure
f = plt.figure()

# create subplots
subplot_dim = (3, 3)
sub11 = plt.subplot2grid(subplot_dim, (0, 0))
sub21 = plt.subplot2grid(subplot_dim, (1, 0), sharex=sub11, sharey=sub11)
sub31 = plt.subplot2grid(subplot_dim, (2, 0), sharex=sub11, sharey=sub11)
sub12 = plt.subplot2grid(subplot_dim, (0, 1))
sub22 = plt.subplot2grid(subplot_dim, (1, 1), sharex=sub12, sharey=sub12)
sub32 = plt.subplot2grid(subplot_dim, (2, 1), sharex=sub12, sharey=sub12)
sub13 = plt.subplot2grid(subplot_dim, (0, 2))
sub23 = plt.subplot2grid(subplot_dim, (1, 2), sharex=sub13, sharey=sub13)
sub33 = plt.subplot2grid(subplot_dim, (2, 2), sharex=sub13, sharey=sub13)

sub11.annotate(r'$\mathbf{q}$ vs $\theta$', xy=(0.5, 1.1), xytext=(0, 0), xycoords='axes fraction', textcoords='offset points', size='large', ha='center', va='baseline')
sub12.annotate(r'$\mathbf{q}$ vs $\phi$', xy=(0.5, 1.1), xytext=(0, 0), xycoords='axes fraction', textcoords='offset points', size='large', ha='center', va='baseline')
sub13.annotate(r'$\theta$ vs $\phi$', xy=(0.5, 1.1), xytext=(0, 0), xycoords='axes fraction', textcoords='offset points', size='large', ha='center', va='baseline')

sub11.annotate('Exact', xy=(0, 0.5), xytext=(-sub11.yaxis.labelpad - 1.1, 0), xycoords=sub11.yaxis.label, textcoords='offset points', size='large', ha='right', va='center', rotation=90)
sub21.annotate('Naive \n n = ' + str(nNaive), xy=(0, 0.5), xytext=(-sub21.yaxis.labelpad - 7, 0), xycoords=sub21.yaxis.label, textcoords='offset points', size='large', ha='center', va='center', rotation=90)
sub31.annotate('Sophisticated \n n = ' + str(nSoph), xy=(0, 0.5), xytext=(-sub31.yaxis.labelpad - 7, 0), xycoords=sub31.yaxis.label, textcoords='offset points', size='large', ha='center', va='center', rotation=90)

plt.suptitle(r'$\mathbf{P}_0 = $' + str(P))


n = 40
qMax = 10

qs = np.linspace(0, qMax, n)
thetas = np.linspace(0, np.pi, n)
phis = np.linspace(0, 2*np.pi, n)
gridQTheta = np.zeros([n, n])
gridQPhi = np.zeros([n, n])
gridThetaPhi = np.zeros([n, n])

feynmanDiagram = FeynmanDiagram(t, P)
v1 = feynmanDiagram.insertVertex(0, dt)
v2 = feynmanDiagram.insertVertex(1, dt)
feynmanDiagram.addInternalPhonon(v1, v2, np.array([1, 1, 1]), 1, 1)
for i, q in enumerate(qs):
  for j, theta in enumerate(thetas):
    for k, phi in enumerate(phis):

      Q = np.array([
        q*np.sin(theta)*np.cos(phi),
        q*np.sin(theta)*np.sin(phi),
        q*np.cos(theta)
      ]);

      feynmanDiagram.setInternalPhononMomentum(feynmanDiagram.Ds[0], Q)
      feynmanDiagram.setInternalPhononMomentumAngle(feynmanDiagram.Ds[0], theta, phi)

      diag = feynmanDiagram()

      gridQTheta[j][i] += diag/n
      gridQPhi[k][i] += diag/n
      gridThetaPhi[k][j] += diag/n


sub11.imshow(gridQTheta, origin='lower', extent=(0, qMax, 0, np.pi))
sub12.imshow(gridQPhi, origin='lower', extent=(0, qMax, 0, 2*np.pi))
sub13.imshow(gridThetaPhi, origin='lower', extent=(0, np.pi, 0, 2*np.pi))



feynmanDiagram = FeynmanDiagram(t, P)
v1 = feynmanDiagram.insertVertex(0, dt)
v2 = feynmanDiagram.insertVertex(1, dt)
feynmanDiagram.addInternalPhonon(v1, v2, np.array([1, 1, 1]), 1, 1)
Da, Ma, Ta, Pa = DMC(feynmanDiagram, nSoph, True)

feynmanDiagram = FeynmanDiagram(t, P)
v1 = feynmanDiagram.insertVertex(0, dt)
v2 = feynmanDiagram.insertVertex(1, dt)
feynmanDiagram.addInternalPhonon(v1, v2, np.array([1, 1, 1]), 1, 1)
Db, Mb, Tb, Pb = DMC(feynmanDiagram, nNaive, False)


sub21.hist2d(Mb, Tb, bins=40, normed=True)
sub22.hist2d(Mb, Pb, bins=40, normed=True)
sub23.hist2d(Tb, Pb, bins=40, normed=True)

sub31.hist2d(Ma, Ta, bins=40, normed=True)
sub32.hist2d(Ma, Pa, bins=40, normed=True)
sub33.hist2d(Ta, Pa, bins=40, normed=True)






















plt.savefig('plot.pdf')

plt.show()