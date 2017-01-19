import numpy as np
from math import factorial

def raiseOrder(self):
  # save old diagram value
  diagOld = self.FD()
  # save diagarm order
  n = len(self.FD.Ds)

  # select two electron propagator on random
  i1, i2 = tuple(np.random.randint(0, len(self.FD.Gs), 2))
  if i1 > i2:
    i1, i2 = i2, i1

  ##
  ## temp
  ##
  # i1, i2 = 0, 0

  g1 = self.FD.Gs[i1]

  # uniformly select a time to split the electron line
  t1 = np.random.uniform(g1.start.position, g1.end.position)
  dt1 = t1 - g1.start.position
  wInvt1 = g1.end.position - g1.start.position

  # split electron
  v1 = self.FD.insertVertex(i1, dt1)

  # after inserting the electron propagator the list got longer
  i2 += 1
  g2 = self.FD.Gs[i2]

  # uniformly select a time to split the electron line
  t2 = np.random.uniform(g2.start.position, g2.end.position)
  dt2 = t2 - g2.start.position
  wInvt2 = g2.end.position - g2.start.position

  # split electron
  v2 = self.FD.insertVertex(i2, dt2)

  # add phonon
  self.FD.addInternalPhonon(v1, v2, np.random.rand(3), 1, 1)

  # take into permutations of i1 and i2 leading to the same diagram
  perm = 1
  if i1 + 1 != i2:
    perm = np.random.randint(0, 2)

  r = factorial(2*n + 1)/(n + 1)


  R = perm * r

  return R


# def raiseOrder(self):
#   # save old diagram value
#   diagOld = self.FD()
#   # save diagarm order
#   n = len(self.FD.Ds)

#   # select two electron propagator on random
#   i1, i2 = tuple(np.random.randint(0, len(self.FD.Gs), 2))
#   if i1 > i2:
#     i1, i2 = i2, i1

#   ##
#   ## temp
#   ##
#   # i1, i2 = 0, 0

#   g1 = self.FD.Gs[i1]

#   # uniformly select a time to split the electron line
#   t1 = np.random.uniform(g1.start.position, g1.end.position)
#   dt1 = t1 - g1.start.position
#   wInvt1 = g1.end.position - g1.start.position

#   # split electron
#   v1 = self.FD.insertVertex(i1, dt1)

#   # after inserting the electron propagator the list got longer
#   i2 += 1
#   g2 = self.FD.Gs[i2]

#   # uniformly select a time to split the electron line
#   t2 = np.random.uniform(g2.start.position, g2.end.position)
#   dt2 = t2 - g2.start.position
#   wInvt2 = g2.end.position - g2.start.position

#   # split electron
#   v2 = self.FD.insertVertex(i2, dt2)


#   # sample phonon momentum
#   std = (t2 - t1)**-0.5;
#   q = np.abs(np.random.normal(0, std))
#   theta = np.random.uniform(0, np.pi)
#   phi = np.random.uniform(0, 2*np.pi)

#   Q = np.array([
#     q*np.sin(theta)*np.cos(phi),
#     q*np.sin(theta)*np.sin(phi),
#     q*np.cos(theta)
#   ])
#   # wInvQ = 2*np.pi**2 * (np.pi/2*std**2)**0.5*np.exp(0.5*(q/std)**2)

#   # add phonon
#   self.FD.addInternalPhonon(v1, v2, Q, theta, phi) 

#   # get current diagram value
#   diag = self.FD()

#   # ration
#   r = (n + 1)/( (2*n + 1)**2 )

#   # R = r*diag/diagOld * wInvt1 * wInvt * wInvQ
#   R = r

#   return R




# def raiseOrder(self):
#   # save old diagram value
#   diagOld = self.FD()
#   # save diagarm order
#   n = len(self.FD.Ds)

#   # select electron propagator on random
#   gIndex = np.random.randint(0, len(self.FD.Gs))
#   g = self.FD.Gs[gIndex]

#   # uniformly select a time to split the electron line
#   t1 = np.random.uniform(g.start.position, g.end.position)
#   dt1 = t1 - g.start.position
#   wInvt1 = g.end.position - g.start.position

#   # split electron
#   v1 = self.FD.insertVertex(gIndex, dt1)

#   # sample phonon length "t" using an exponential distribution
#   dt = self.FD.time - t1
#   r = np.random.rand()
#   t = -np.log(1 - r*(1 - np.exp(-dt)))
#   wInvt = (1 - np.exp(-dt))*np.exp(t)

#   # look for corresponding electron propagator to split
#   i = gIndex + 1
#   g = self.FD.Gs[i]
#   while g.end.position < t1 + t and g.end.G[1]:
#     i += 1
#     g = g.end.G[1]
#   v2 = self.FD.insertVertex(i, t1 + t - g.start.position)

#   ##
#   ## chose on random propagator after v1 to split
#   # ##
#   # if True:
#   #   gIndex = np.random.randint(gIndex + 1, len(self.FD.Gs))
#   #   g = self.FD.Gs[gIndex]

#   #   rtemp = len(self.FD.Gs) - 1 - (gIndex + 1)

#   #   t2 = np.random.uniform(g.start.position, g.end.position)
#   #   dt2 = t2 - g.start.position

#   #   # split electron
#   #   v2 = self.FD.insertVertex(gIndex, dt2)


#   # sample phonon momentum
#   std = t**-0.5;
#   q = np.abs(np.random.normal(0, std))
#   theta = np.random.uniform(0, np.pi)
#   phi = np.random.uniform(0, 2*np.pi)

#   Q = np.array([
#     q*np.sin(theta)*np.cos(phi),
#     q*np.sin(theta)*np.sin(phi),
#     q*np.cos(theta)
#   ])
#   # wInvQ = 2*np.pi**2 * (np.pi/2*std**2)**0.5*np.exp(0.5*(q/std)**2)

#   # add phonon
#   self.FD.addInternalPhonon(v1, v2, Q, theta, phi) 

#   # get current diagram value
#   diag = self.FD()

#   # ration
#   r = (n + 1)/( (2*n + 1) )

#   # R = r*diag/diagOld * wInvt1 * wInvt * wInvQ
#   # R = r

#   # R = r * wInvt1 * wInvt

#   R = 1

#   return R