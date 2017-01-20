import numpy as np

def swapPhononConnections(self):
  # pick an internal electron propagator on random
  g = np.random.choice(self.FD.Gs[1:-1])

  v1 = g.start
  v2 = g.end

  c1 = 1 if v1.D[1] else -1
  c2 = 1 if v2.D[1] else -1

  d1 = v1.D[0] or v1.D[1]
  d2 = v2.D[0] or v2.D[1]

  t = v2.position - v1.position
  Eafter = 0.5*np.linalg.norm((g.momentum + c1*d1.momentum - c2*d2.momentum))**2
  Ebefore = 0.5*np.linalg.norm(g.momentum)**2

  # if both vertices belong to the same phonon
  # same amount of phonon is still present above after the swap in this case
  if d1 == d2:
    dw = 0
  else:
    dw = c2 - c1

  exponent = -t*(Eafter - Ebefore + dw)
  if exponent > 700:
    # to prevent overflow. could have set the treshhold to 0 since R >= 1 all do the same thing
    exponent = 700
  
  R = np.exp(exponent)

  if self.debug:
    diagOld = self.FD()

  self.FD.swapPhononEnds(v1, v2)

  if self.debug:
    diag = self.FD()

    if exponent == 700:
      print('--------------------------------------------------------------------')
      print('Connection Swap overflow', R*diagOld/diag, 'diag=', diag, 'diagOld=', diagOld)
      print('--------------------------------------------------------------------')
    else:
      print('Connection Swap', R*diagOld/diag)

  R = 1

  return R