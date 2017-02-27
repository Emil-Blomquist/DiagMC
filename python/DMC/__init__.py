import sys
import numpy as np

sys.path.append('../')
from FeynmanDiagram import FeynmanDiagram

class DiagrammaticMonteCarlo(object):

  sys.path.append('/updates')

  from .updates.changeInternalPhononMomentumDirection import changeInternalPhononMomentumDirection
  from .updates.changeInternalPhononMomentumMagnitude import changeInternalPhononMomentumMagnitude
  from .updates.swapPhononConnections import swapPhononConnections
  from .updates.shiftVertexPosition import shiftVertexPosition
  from .updates.raiseOrder import raiseOrder
  from .updates.lowerOrder import lowerOrder

  def __init__(self, params, debug = False):
    self.t = params['t']
    self.mu = params['mu']
    self.a = params['a']
    self.P = params['P']
    self.maxOrder = params['maxOrder']
    self.minOrder = params['minOrder']

    # debug mode
    self.debug = debug

    ##  
    ## in the future we should insert "params" into FeynmanDiagram
    ##

    # initiate Feynman diagram
    self.FD = FeynmanDiagram(self.t, self.P)

    ##
    ## in future we should let the actual algorithm build the diagram
    ##
    dt = 0.1;

    v1 = self.FD.insertVertex(0, dt)
    v2 = self.FD.insertVertex(1, dt)
    v3 = self.FD.insertVertex(2, dt)
    v4 = self.FD.insertVertex(3, dt)
    # v5 = self.FD.insertVertex(4, dt)
    # v6 = self.FD.insertVertex(5, dt)

    self.FD.addInternalPhonon(v1, v3, np.array([1, 2, 1]), 1, 2) 
    self.FD.addInternalPhonon(v2, v4, np.array([3, 4, 5]), 1, 2)
    # self.FD.addInternalPhonon(v4, v6, np.random.rand(3), 1, 1)

  def run(self, N):
    ##
    ## here the actual DMC takes place
    ##

    updates = [
      self.changeInternalPhononMomentumDirection,
      # self.changeInternalPhononMomentumMagnitude,
      # self.shiftVertexPosition,
      # self.swapPhononConnections,
      # self.raiseOrder,
      # self.lowerOrder
    ]
    

    # to reach some sense of randomness
    # for i in range(0, 10):
    #   for update in updates:
    #     if len(self.FD.Ds) == 0 and update != self.raiseOrder:
    #       r = 1
    #     elif len(self.FD.Ds) == self.minOrder and (update == self.lowerOrder):
    #       r = 1
    #     elif len(self.FD.Ds) == self.maxOrder and (update == self.raiseOrder):
    #       r = 1
    #     else:
    #       update()
   

    diagBins = {}
    orderBins = {}
    for i in range(0, N):
      self.FD.save()

      update = np.random.choice(updates)

      if len(self.FD.Ds) == 0 and update != self.raiseOrder:
        r = 1
      elif len(self.FD.Ds) == self.minOrder and (update == self.lowerOrder):
        r = 1
      elif len(self.FD.Ds) == self.maxOrder and (update == self.raiseOrder):
        r = 1
      else:
        r = update()

      # r = update()

      # if np.random.rand() > min(1, r):
      #   # reject step
      #   self.FD.revert()

      self.FD.plot()

      ##
      ## temp
      ##
      # else: 
      #   struct = self.FD.structure()
      #   if struct in diagBins:
      #     diagBins[struct] += 1
      #   else:
      #     diagBins[struct] = 1

      #   self.FD.revert()


      struct = self.FD.structure()
      if struct in diagBins:
        diagBins[struct] += 1
      else:
        diagBins[struct] = 1

      order = len(self.FD.Ds)
      if order in orderBins:
        orderBins[order] += 1
      else:
        orderBins[order] = 1

    for diag in diagBins:
      diagBins[diag] = diagBins[diag]/N

    for order in orderBins:
      orderBins[order] = orderBins[order]/N

    return diagBins, orderBins