import numpy as np
import matplotlib.pyplot as plt


file1 = 'p=0.000 a=2.000 mu=-2.200 N=20000000 time=2017-02-23 07:05:08'
file2 = 'temp'


T1 = [];
G1 = [];
with open('data/' + file1 + '.txt') as fp:
  for i, line in enumerate(fp):

    # skip first line
    if i == 0:
      print(line[0:-1])
      continue

    l = line.split('\n')[0]
    l = l.split(' ')

    if i == 1:
      Gs1 = np.array(l[1:-1])
    else:
      Gs1 = np.vstack([Gs1, l[1:-1]])

    T1.append(l[0])
    G1.append(l[-1])


T2 = [];
G2 = [];
with open('data/' + file2 + '.txt') as fp:
  for i, line in enumerate(fp):

    # skip first line
    if i == 0:
      continue

    l = line.split('\n')[0]
    l = l.split(' ')

    if i == 1:
      Gs2 = np.array(l[1:-1])
    else:
      Gs2 = np.vstack([Gs2, l[1:-1]])

    T2.append(l[0])
    G2.append(l[-1])





if len(T1) != len(T2):
  print('Not the same amount of data points')
else:

  for i in range(0, len(T1)):
    dt = abs(float(T1[i]) - float(T2[i]))
    if dt > 10**-5:
      print('ERROR: times do not match at row ', i, ' which corresponds to t≈', T1[i])
    else:
      string = "{:.7f}".format(0.5*(float(T1[i]) + float(T2[i]))) + ' '
      for j in range(0, Gs1.shape[1]):
        string += "{:.7f}".format(0.5*(float(Gs1[i, j]) + float(Gs2[i, j]))) + ' '
      string += "{:.7f}".format(0.5*(float(G1[i]) + float(G2[i])))

    print(string)

# plot G
plt.plot(T1, G1, lw=2, color="brown", label=r'$\sum G^i$')
plt.plot(T2, G2, lw=2, color="red", label=r'$\sum G^i$')
plt.show()