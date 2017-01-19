import numpy as np

def calculateQ(P0, q, theta, phi):
  p0 = np.linalg.norm(P0)

  # 10^-10 to handle rounding errors
  if p0 > 10**-10:
    Ep = P0/p0
  else:
    # if P0 ≈ (0,0,0) we use that theta is the angle against the z-axis
    Ep = np.array([0, 0, 1])

  # generate a temporary vector used to obtain two other vectors in order to span the rest of R^3
  tempVector = Ep + np.array([1, 0, 0])
  Eo1 = np.cross(Ep, tempVector)
  if np.linalg.norm(Eo1) < 10**-14:
    tempVector = Ep + np.array([0, 1, 0])
    Eo1 = np.cross(Ep, tempVector)

  # two orthogonal vectors spanning the plane normal to Ep
  Eo1 = Eo1/np.linalg.norm(Eo1)
  Eo2 = np.cross(Ep, Eo1)

  # generating parallel and orthogonal component of Q
  Qp = Ep * q*np.cos(theta)
  Qo = (Eo1*np.cos(phi) + Eo2*np.sin(phi)) * q*np.sin(theta)

  return Qp + Qo