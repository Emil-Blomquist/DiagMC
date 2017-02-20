import numpy as np
import matplotlib.pyplot as plt

T = [];
G = [];

if False:
  file = "data-0-0-0.txt"
  P = "(0, 0, 0)"
else:
  file = "data-1-0-0.txt"
  P = "(1, 0, 0)"


with open(file) as fp:
  for i, line in enumerate(fp):

    l = line.split('\n')[0]
    l = l.split(' ')

    if i == 0:
      Gs = np.array(l[1:-1])
    else:
      Gs = np.vstack([Gs, l[1:-1]])

    # print(l)

    T.append(l[0])
    G.append(l[-1])

for i, val in enumerate(Gs[0, :]):
  plt.plot(T, Gs[:, i], label=r'$G^' +  str(i) + '$')

plt.plot(T, G, lw=2, color="brown", label=r'$\sum G^i$')

plt.xlim(0, 5)
plt.xlabel(r'$\tau$', fontsize=20)

titleText = '$\mathbf{p} = ' + P + ', \, \alpha = 2, \, \mu = -2.2$'

plt.title(r'$\mathbf{{p}} = {0}, \, \alpha = 2, \, \mu = -2.2$'.format(P))
plt.legend(loc=1)

plt.savefig('plot.pdf')
plt.show()