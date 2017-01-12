def structure(self):
  Ds = dict()
  name = ''

  for g in self.Gs[0:-1]:
    v = g.end

    if v.D[1]:
      # outgoing
      name += str(len(Ds) + 1)
      Ds[id(v.D[1])] = len(Ds) + 1
    else:
      # incoming
      i = Ds[id(v.D[0])]
      name += str(i)

  return int(name)