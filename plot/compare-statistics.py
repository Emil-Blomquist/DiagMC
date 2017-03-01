import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# fileName = 'p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 *'
# paths = glob.glob('data/determine lambda exp/' + fileName + '.txt', recursive=True)


paths = [
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=0 param=0.100.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=1 param=0.270.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=2 param=0.439.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=3 param=0.609.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=4 param=0.778.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=5 param=0.948.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=6 param=1.117.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=7 param=1.287.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=8 param=1.457.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=9 param=1.626.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=10 param=1.796.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=11 param=1.965.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=12 param=2.135.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=13 param=2.304.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=14 param=2.474.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=15 param=2.643.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=16 param=2.813.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=17 param=2.983.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=18 param=3.152.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=19 param=3.322.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=20 param=3.491.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=21 param=3.661.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=22 param=3.830.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 12:35:16 id=23 param=4.000.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=0 param=4.000.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=1 param=4.170.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=2 param=4.339.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=3 param=4.509.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=4 param=4.678.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=5 param=4.848.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=6 param=5.017.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=7 param=5.187.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=8 param=5.357.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=9 param=5.526.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=10 param=5.696.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=11 param=5.865.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=12 param=6.035.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=13 param=6.204.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=14 param=6.374.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=15 param=6.543.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=16 param=6.713.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=17 param=6.883.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=18 param=7.052.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=19 param=7.222.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=20 param=7.391.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=21 param=7.561.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=22 param=7.730.txt',
  'data/determine lambda exp/p=0.000 a=2.000 mu=-2.200 N=1000000 time=2017-02-28 16:56:07 id=23 param=7.900.txt'
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

    plt.plot(config['param'], np.std(G), '.', color=colors[j%len(colors)])

plt.savefig('plots/' + 'many' + '.pdf')
plt.show()