from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def momentaExtremum(Gs):
  P = []

  # not very optimised algorithm since P will contain multiples
  for i in range(0, len(Gs)):
    v = Gs[i].end

    if v.G[0]:
      P.append(norm(v.G[0].momentum))
    if v.G[1]:
      P.append(norm(v.G[1].momentum))

    if v.D[0]:
      P.append(norm(v.D[0].momentum))
    elif v.D[0]:
      P.append(norm(v.D[0].momentum))

  return np.min(P), np.max(P)

def roundNumber(num):
  return str(np.round(num*1000)/1000)

def norm(num):
  if np.size(num) == 1:
    return num
  else:
    return np.linalg.norm(num)

def colorCode(pMin, pMax, P):
  if pMin == pMax:
    c = [0, 0, 1]
  else:
    r = (norm(P) - pMin)/(pMax - pMin)
    c = [r, 0, 1 - r]

  return c

def plot(self):
  pMin, pMax = momentaExtremum(self.Gs)

  # create a list of times
  # is also used to identify vertices for D's
  times = [0]
  for g in self.Gs:
    times.append(g.end.position)

  # figure settings. Use ax for plot
  # plt.ion()
  fig = plt.figure(facecolor='white', tight_layout=True, figsize=(1.5*len(self.Gs), 1.5*0.5*len(self.Gs)))
  plt.axis('off')
  ax = plt.gca()
  ax.set_aspect('equal')
  ax.margins(0.05)

  ##
  ## How to include external phonons?
  ##


  ##
  ## Check for validity of diagram
  ##
  for i, g in enumerate(self.Gs[1:], 1):
    v = g.start
    invalid = False

    # Three lines
    if not (v.G[0] and v.G[1] and bool(v.D[0]) ^ bool(v.D[1])):
      invalid = True
      print(i, 'DIAGRAM ERROR: Not exactly three lines connected to a vertex')
      for g in v.G + v.D:
        startTime = str(roundNumber(g.start.position)) if g.start else '?'
        endTime = str(roundNumber(g.end.position)) if g.start else '?'
        print('\t' + g.type + ': ' + startTime + ' -> ' + endTime, g.momentum)

    # Conservation of momentum
    dP = v.G[1].momentum - v.G[0].momentum
    if v.D[1]:
      dP += v.D[1].momentum
    elif v.D[0]:
      dP -= v.D[0].momentum

    if norm(dP) > 10**-10:
      invalid = True
      print(i, 'DIAGRAM ERROR: Momentum nonconservation over a vertex')
      print('\tdp =', norm(dP))

    # Propagators connected to two vertices
    ##
    ## should probably do this for G and D and not in vertices validity checker
    ## however sine it has a loose end it will not be drawn...
    ## perhaps draw loose ends to -1?
    # if 3 != len([g for g in v.G + v.D if g is not None]):
    #   print(i, 'DIAGRAM ERROR: Propagator with loose end found')
    #   invalid = True

    # increasing times for vertices
    if not v.G[0].start.position < v.position < v.G[1].end.position:
      print(i, 'DIAGRAM ERROR: Vertices not chronologically ordered')
      invalid = True

    # mark erroneous vertex
    if invalid:
      r = 0.2
      a = np.linspace(0, 2*np.pi, 100)
      x = i + r*np.cos(a)
      y = r*np.sin(a)

      plt.fill(x, y, zorder=1, c=np.array([42, 252, 92])/255)

  ##
  ## check the phonon propagators
  ##
  for d in self.Ds:
    invalid = False

    # propagating backwards
    if d.start.position > d.end.position:
      print('DIAGRAM ERROR: Phonon propagating backwards')
      invalid = True

    # mark invalid phonon propagators
    if invalid:
      i = times.index(d.start.position)
      j = times.index(d.end.position)

      # plot dotted arc
      a = np.linspace(0, np.pi, 100)
      x = 0.5*(i + j) + 0.5*(j - i)*np.cos(a)
      y = 0.5*(j - i)*np.sin(a)
      ax.plot(x, y, lw=15, zorder=1, c=np.array([42, 252, 92])/255)

  ##
  ## plot diagram
  ##

  # plot vertices
  for i, time in enumerate(times):
    # color
    color = 'black'
    if time == 0 or time == self.time:
      color = 'white'

    # plot vertex dot
    ax.scatter(i, 0, s=40, c=color, lw=1.5, zorder=2)
    # plot vertex time
    ax.text(i, -0.15, roundNumber(time), horizontalalignment='center')

  # plot Gs
  for i, g in enumerate(self.Gs):
    # plot solid line
    ax.plot([i, i + 1], [0, 0], lw=2, zorder=1, c=colorCode(pMin, pMax, g.momentum))
    # plot G momentum
    ax.text(i + 0.5, 0.05, roundNumber(norm(g.momentum)), horizontalalignment='center')

  # plot Ds
  for d in self.Ds:
    i = times.index(d.start.position)
    j = times.index(d.end.position)

    # plot dotted arc
    a = np.linspace(0, np.pi, 100)
    x = 0.5*(i + j) + 0.5*(j - i)*np.cos(a)
    y = 0.5*(j - i)*np.sin(a)
    ax.plot(x, y, '--', lw=2, zorder=1, c=colorCode(pMin, pMax, d.momentum))
    # plot D momentum
    ax.text(0.5*(i + j), 0.5*(j - i) + 0.05, roundNumber(norm(d.momentum)), horizontalalignment='center',
      bbox={'facecolor':'white', 'pad':0, 'lw':0}, zorder=2)


  plt.show()
  plt.close()