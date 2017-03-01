import numpy as np
import matplotlib.pyplot as plt

fileName = 'G7 p=0 a=2'
T = [];
G = [];

with open('data/' + fileName + '.txt') as fp:
  for i, line in enumerate(fp):

    # fetch data from first line
    if i == 0:
      data = line[9: -10].split(' ')
      config = {}
      for d in data:
        separated = d.split('=', 1)
        if (len(separated) == 1):
          config[key] += ' ' + separated[0]
        else:
          key = separated[0]
          config[key] = separated[-1]
    else:

      l = line.split('\n')[0]
      l = l.split(' ')

      if i == 1:
        Gs = np.array(l[1:-1])
      else:
        Gs = np.vstack([Gs, l[1:-1]])

      T.append(l[0])
      G.append(l[-1])

  T = [float(i) for i in T]
  G = [float(i) for i in G]

for i, val in enumerate(Gs[0, :]):
  plt.plot(T, Gs[:, i], label=r'$G^' +  str(i) + '$')

print(G)
print(Gs[:, 0])
print(Gs[:, 1])

plt.plot(T, G, lw=2, color="brown", label=r'$\sum G^i$')

plt.xlim(0, 5)
plt.xlabel(r'$\tau$', fontsize=20)
plt.ylabel(config['time'])

plt.title(r'$p = {0}, \, \alpha = {1}, \, \mu = {2}, \, N = {3}$'.format(config['p'], config['a'], config['mu'], config['N']))
plt.legend(loc=1)

plt.savefig('plots/' + fileName + '.pdf')
plt.show()