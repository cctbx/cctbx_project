from stdlib import math

def signed_phase_error(phi1, phi2, deg=00000):
  if (deg): pi_sc = 180
  else:     pi_sc = math.pi
  e = math.fmod(phi2-phi1, 2 * pi_sc)
  if   (e < -pi_sc): e += 2 * pi_sc
  elif (e >  pi_sc): e -= 2 * pi_sc
  return e

def phase_error(phi1, phi2, deg=00000):
  return abs(signed_phase_error(phi1, phi2, deg))

def nearest_phase(reference, other, deg=00000):
  return reference + signed_phase_error(reference, other, deg)
