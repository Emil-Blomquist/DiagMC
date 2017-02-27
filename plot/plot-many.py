import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

T = [];
G = [];

fileName = 'p=0.000 a=2.000 mu=-2.200 N=40000000 *'
# fileName = 'temp'

paths = glob.glob('data/' + fileName + '.txt', recursive=True)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for path in paths:
  print(path)


  T = [];
  G = [];

  with open(path) as fp:
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

    colors = ['cyan', 'purple', 'orange', 'red', 'magenta', 'blue', 'black', 'gray']

    # for i in range(10, len(Gs[0, :])):
      # ax.plot(T, np.full(len(T), float(config['p'])), Gs[:, i].astype(np.float), color=colors[i])
    
    ax.plot(T, np.full(len(T), float(config['p'])), G, color='brown')


plt.savefig('plots/' + 'many' + '.pdf')
plt.show()