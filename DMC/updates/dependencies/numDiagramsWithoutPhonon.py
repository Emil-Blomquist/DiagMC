
def numDiagramsWithoutPhonon(FD):
  
  names = dict()

  if len(FD.Ds) == 1:
    return 1

  for d in FD.Ds:

    Ds = dict()
    name = ''
    for g in FD.Gs[0:-1]:
      v = g.end

      if v.D[1] and v.D[1] != d:
        # outgoing
        name += str(len(Ds) + 1)
        Ds[id(v.D[1])] = len(Ds) + 1
      elif v.D[0] and v.D[0] != d:
        # incoming
        i = Ds[id(v.D[0])]
        name += str(i)

    names[name] = True

  return len(names)
