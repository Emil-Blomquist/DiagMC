import numpy as np
import matplotlib.pyplot as plt


a = 100

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

