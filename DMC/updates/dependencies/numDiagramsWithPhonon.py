
def numDiagramsWithPhonon(FD):
  
  names = dict()

  if len(FD.Ds) == 0:
    return 1

  for g1 in FD.Gs[0:-1]:
    g2 = g1.end.G[1]

    Ds = dict()
    name = ''
    for g in FD.Gs[0:-1]:
      v = g.end

      if g == g1:
        # outgoing
        name += str(len(Ds) + 1)
        Ds['artificial'] = len(Ds) + 1

      if v.D[1]:
        # outgoing
        name += str(len(Ds) + 1)
        Ds[id(v.D[1])] = len(Ds) + 1
      elif v.D[0]:
        # incoming
        i = Ds[id(v.D[0])]
        name += str(i)

      if g.end.G[1] == g2:
        # incoming
        i = Ds['artificial']
        name += str(i)

    names[name] = True

  return len(names)
