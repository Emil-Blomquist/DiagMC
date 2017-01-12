from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

from FeynmanDiagram import FeynmanDiagram

##
##
## If we use the proper probability distribution when sampling
## internal variables we wont need to save the invProb in the
## feynmanDiagram object since all of that will be taken care of here
##
##

##
## shift phonon position
##
def shiftPhononPosition(fd):
  # pick an internal phonon propagator on random
  d = np.random.choice(fd.Ds)
  v1 = d.start
  v2 = d.end

  # first we shift the front verex
  t11 = v1.G[0].start.position
  t12 = v1.G[1].end.position
  t1 = np.random.uniform(t11, t12)

  invProb1 = t12 - t11
  fd.setVertexPosition(v1, t1, invProb1)

  # then the back one
  t21 = v2.G[0].start.position
  t22 = v2.G[1].end.position
  t2 = t21 - np.log(1 - np.random.rand()*(1 - np.exp(-t22 + t21)))

  invProb2 = (1 - np.exp(t21 - t22))/np.exp(t21 - t2)
  fd.setVertexPosition(v2, t2, invProb2)

##
## swap phonon ends
##
def swapPhononEnds(fd):
  # pick an internal electron propagator on random
  g = np.random.choice(fd.Gs[1:-1])

  v1 = g.start
  v2 = g.end

  fd.swapPhononEnds(v1, v2)


##
## change phonon momentum
##
def changePhononMomentum(fd):
  # pick an internal phonon propagator on random
  d = np.random.choice(fd.Ds)

  standDev = (d.end.position - d.start.position)**-0.5;
  q = np.abs(np.random.normal(0, standDev));

  theta = np.random.uniform(0, np.pi);
  phi = np.random.uniform(0, 2*np.pi);

  Q = np.array([
    q*np.sin(theta)*np.cos(phi),
    q*np.sin(theta)*np.sin(phi),
    q*np.cos(theta)
  ]);

  invProb = np.pi**2 * (2*np.pi*standDev**2)**0.5*np.exp(0.5*(q/standDev)**2)

  fd.setInternalPhononMomentum(d, Q, np.sin(theta), invProb)

def DMC(t, N, n):

  # create second order diagram
  feynmanDiagram = FeynmanDiagram(t, np.array([0, 0, 0]))

  v1 = feynmanDiagram.insertVertex(0, 0.001, 1)
  v2 = feynmanDiagram.insertVertex(1, 0.001, 1)
  v3 = feynmanDiagram.insertVertex(2, 0.001, 1)
  v4 = feynmanDiagram.insertVertex(3, 0.001, 1)
  feynmanDiagram.addInternalPhonon(v1, v4, 2*np.array([1, 1, 1]), 1, 1)
  feynmanDiagram.addInternalPhonon(v2, v3, 3*np.array([1, 1, 1]), 1, 1)

  print(n)

  # to reach some sense of randomness
  if n > 0:
    for i in range(0, n):
      shiftPhononPosition(feynmanDiagram)
      swapPhononEnds(feynmanDiagram)
      changePhononMomentum(feynmanDiagram)

  wInv = 1
  bins = {}
  for i in range(0, N):
    feynmanDiagram.save()

    (fdOld, invProbOld) = feynmanDiagram()

    update = np.random.choice([shiftPhononPosition, swapPhononEnds, changePhononMomentum])
    update(feynmanDiagram)

    (fd, invProb) = feynmanDiagram()

    r = fd*invProb/(fdOld*invProbOld)
    a = min(1, r)

    if np.random.rand() > a:
      # reject step
      feynmanDiagram.revert()


    struct = feynmanDiagram.structure()
    if struct in bins:
      bins[struct] += 1
    else:
      bins[struct] = 1

  for key in bins:
    bins[key] = round(1000*bins[key]/N)/1000

  return bins



n = [0] * 10
n = n + ([10] * 10)
n = n + ([100] * 10)
n = n + ([1000] * 10)

plt.ion()

N = 10000
t1122 = []
t1212 = []
t1221 = []
for i in n:
  bins = DMC(5, N, i)

  t1122.append(bins[1122])
  t1212.append(bins[1212])
  t1221.append(bins[1221])

  print(i, bins[1122], bins[1212], bins[1221])

  plt.plot(t1122, 'b.-')
  plt.plot(t1212, 'r.-')
  plt.plot(t1221, 'g.-')
  plt.pause(0.05)

plt.ioff()
plt.show()