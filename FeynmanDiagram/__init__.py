import sys
import numpy as np

sys.path.append('../')
from Vertex import Vertex
from Propagator import G, D


class FeynmanDiagram(object):

  from .plot import plot
  from .phononMethods import addInternalPhonon, removeInternalPhonon, setInternalPhononMomentum, setInternalPhononMomentumAngle, swapPhononEnds
  from .vertexMethods import insertVertex, removeVertex, setVertexPosition
  from .structure import structure

  def __init__(self, time, momentum):
    self.time = time
    self.momentum = momentum

    # list of electron propagator in chronological order
    self.Gs = []
    # list of internal phonons propagators in the order they were created
    self.Ds = []

    # saved variables
    self.timeSaved = None
    self.momentumSaved = None
    self.GsSaved = None
    self.DsSaved = None

    # start and end points
    self.start = Vertex(0)
    self.end = Vertex(time)

    ##
    ## creates a bare propagator
    ## in future we might want to remake this
    ##
    g = G(momentum)
    g.setStart(self.start)
    g.setEnd(self.end)
    self.Gs.append(g)

  def __call__(self):
    ##
    ## include external phonons
    ##

    greens = 1

    ##
    ## not necesarry
    ##
    # chemical potential
    # mu = -2.2
    # greens *= np.exp(mu*self.time)

    # electrons
    for g in self.Gs:
      greens *= g()

    # phonons
    for d in self.Ds:
      greens *= d()


    return greens

  def save(self):
    self.start.save()
    for g in self.Gs:
      g.save()
      g.end.save()

    for d in self.Ds:
      d.save()

    self.timeSaved = self.time
    self.momentumSaved = self.momentum
    self.GsSaved = self.Gs[:]
    self.DsSaved = self.Ds[:]

  def revert(self):
    self.time = self.timeSaved
    self.momentum = self.momentumSaved
    self.Gs = self.GsSaved
    self.Ds = self.DsSaved

    self.start.revert()
    for g in self.Gs:
      g.revert()
      g.end.revert()

    for d in self.Ds:
      d.revert()
