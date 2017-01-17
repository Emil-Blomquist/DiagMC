import numpy as np
import matplotlib.pyplot as plt
import math

from scipy.special import erfinv, erf

##
## W(theta)
##
if  False:
  a = 1

  Theta = np.linspace(0, np.pi, 1000)
  W = 0.5*a*np.sin(Theta)*np.exp(a*np.cos(Theta))/np.sinh(a)
  plt.plot(Theta, W)


  C = []
  for i in range(0, 10000):
    r = np.random.rand()
    cosTheta = 1 + np.log(1 - r*(1 - np.exp(-2*a)))/a
    theta = np.arccos(cosTheta)

    C.append(theta)

  plt.hist(C, 100, normed=True)
  plt.xlim(0, np.pi)
  plt.show()

##
## W(q)
##
if True:
  a = 0.1
  b = 3

  C = 2*(a/np.pi)**0.5 * 1/(1 - math.erf(-a**0.5*b))

  q = np.linspace(0, 10, 1000)
  W = []
  for i in q:
    W.append(C * math.exp(-a*(i - b)**2))


  # scipy.special.erfinv(y)

  C = []
  for i in range(0, 100000):
    r = np.random.rand()
    qq = b + a**-0.5 * erfinv(r + (r - 1)*erf(a**0.5*b))

    C.append(qq)

  plt.hist(C, 100, normed=True)

  plt.plot(q, W)
  plt.show()

