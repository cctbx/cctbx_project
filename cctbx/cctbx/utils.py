import math

def phase_error(phi1, phi2, deg=False):
  if (deg): pi_sc = 180
  else:     pi_sc = math.pi
  e = math.fmod(phi1-phi2, 2 * pi_sc)
  if   (e < -pi_sc): e += 2 * pi_sc
  elif (e >  pi_sc): e -= 2 * pi_sc
  return abs(e)
