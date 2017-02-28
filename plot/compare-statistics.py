import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fileName = 'p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=* param=*'

# paths = glob.glob('data/' + fileName + '.txt', recursive=True)

paths = [
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=1 param=0.250.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=2 param=0.375.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=3 param=0.500.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=4 param=0.625.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=5 param=0.750.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=6 param=0.875.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=7 param=1.000.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=8 param=1.125.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=9 param=1.250.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=10 param=1.375.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=11 param=1.500.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=12 param=1.625.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=13 param=1.750.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=14 param=1.875.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=15 param=2.000.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=16 param=2.125.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=17 param=2.250.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=18 param=2.375.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=19 param=2.500.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=20 param=2.625.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=21 param=2.750.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=22 param=2.875.txt',
  'data/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-27 23:21:46 id=23 param=3.000.txt',
]

for j, path in enumerate(paths):
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

    # for i in range(0, len(Gs[0, :])):
    #   plt.plot(T, Gs[:, i].astype(np.float), color=colors[j])

    # plt.plot(G, color=colors[j%len(colors)])

    plt.plot(j, np.std(G), '.', color=colors[j%len(colors)])

plt.savefig('plots/' + 'many' + '.pdf')
plt.show()