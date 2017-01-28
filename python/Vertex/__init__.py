from __future__ import division
import sys, traceback
import numpy as np
import matplotlib.pyplot as plt

class Vertex:
  def __init__(self, pos):
    self.position = pos

    self.G = [None, None]
    self.D = [None, None]

    # saved variables
    self.positionSaved = None
    self.GSaved = None
    self.DSaved = None

  def setPropagatorIn(self, prop):
    if prop.type == 'G':
      self.G[0] = prop
    else:
      self.D[0] = prop

  def setPropagatorOut(self, prop):
    if prop.type == 'G':
      self.G[1] = prop
    else:
      self.D[1] = prop

  def setPosition(self, position):
    self.position = position;

  def save(self):
    self.positionSaved = self.position
    self.GSaved = self.G[:]
    self.DSaved = self.D[:]

  def revert(self):
    self.position = self.positionSaved
    self.G = self.GSaved
    self.D = self.DSaved
