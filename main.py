from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

from FeynmanDiagram import FeynmanDiagram



##
## shift vertex position
##
def shiftVertexPosition(fd):
  # pick a vertex on random
  g = np.random.choice(fd.Gs[0:-1])
  v = g.end

  # fetch available time interval
  t1 = g.start.position
  t2 = v.G[1].end.position
  t = np.random.uniform(t1, t2)

  fd.setVertexPosition(v, t)

##
## swap phonon ends
##
# def swapPhononEnd(fd):
  ##
  ## for this process we do NOT have detailed balance.
  ##

  # ##
  # ## we must not allow  vertices to belong to the same propagator (detailed balance broken)
  # ##
  # if len(fd.Ds) < 2:
  #   print('remove swapPhononEnd from possible updates')

  # # pick an internal electron propagator on random
  # Gs = fd.Gs[1:-1]
  # while True:
  #   gIndex = np.random.randint(0, len(Gs))
  #   g = Gs[gIndex]
  #   v1 = g.start
  #   v2 = g.end

  #   d1 = v1.D[0] or v1.D[1]
  #   d2 = v2.D[0] or v2.D[1]

  #   if d1 == d2:
  #     del Gs[gIndex]
  #   else:
  #     break


  # fd.swapPhononEnds(v1, v2)
def swapPhononEnd(fd):
    # pick an internal electron propagator on random
  g = np.random.choice(fd.Gs[1:-1])

  v1 = g.start
  v2 = g.end

  c1 = 1 if v1.D[1] else -1
  c2 = 1 if v2.D[1] else -1

  d1 = v1.D[0] or v1.D[1]
  d2 = v2.D[0] or v2.D[1]

  # if both vertices belong to the same phonon
  if v1.D[1] == v2.D[0]:
    c1 = c2 = 0

  t = v2.position - v1.position
  Eafter = 0.5*(g.momentum + c1*d1.momentum - c2*d2.momentum)**2
  Ebefore = 0.5*g.momentum**2
  R = np.exp(-t*(Eafter - Ebefore -(c1 - c2)))

  fd.swapPhononEnds(v1, v2)

  return R


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

  fd.setInternalPhononMomentum(d, Q, np.sin(theta))









feynmanDiagram = FeynmanDiagram(1, np.array([0, 0, 0]))

v1 = feynmanDiagram.insertVertex(0, 0.1, 1)
v2 = feynmanDiagram.insertVertex(1, 0.1, 1)
v3 = feynmanDiagram.insertVertex(2, 0.1, 1)
v4 = feynmanDiagram.insertVertex(3, 0.1, 1)
feynmanDiagram.addInternalPhonon(v1, v4, 2*np.array([1, 1, 1]), 1, 1)
feynmanDiagram.addInternalPhonon(v2, v3, 3*np.array([1, 1, 1]), 1, 1)


# feynmanDiagram.save()
# feynmanDiagram.plot()
# print(feynmanDiagram())


# v5 = feynmanDiagram.insertVertex(4, 0.1)
# v6 = feynmanDiagram.insertVertex(5, 0.1)
# v7 = feynmanDiagram.insertVertex(6, 0.1)
# v8 = feynmanDiagram.insertVertex(7, 0.1)
# feynmanDiagram.addInternalPhonon(v5, v6, 3*np.array([1, 1, 1]), 1)
# feynmanDiagram.addInternalPhonon(v7, v8, 3*np.array([1, 1, 1]), 1)


bins = {}
for i in range(0, 1000):
  swapPhononEnd(feynmanDiagram)

  # name = feynmanDiagram.structure()
  # if name in bins:
  #   bins[name] += 1
  # else:
  #   bins[name] = 1


  feynmanDiagram.plot()

feynmanDiagram.plot()
print(feynmanDiagram())
feynmanDiagram.revert()
feynmanDiagram.plot()
print(feynmanDiagram())


# print(feynmanDiagram())

print(len(bins), bins)
plt.plot(list(bins.values()))
plt.show()