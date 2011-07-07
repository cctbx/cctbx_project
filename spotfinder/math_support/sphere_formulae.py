import math

"""radius is taken to be d_star, or 1/resolution limit in units of distance"""

def sphere_volume(radius):
  return (4./3.)*math.pi*math.pow(radius,3)

def sph_volume_minus_missing_cone(radius,wavelength):
  R = radius
  K = 1./wavelength
  radical = math.sqrt(K*K - R*R)
  return 2. * math.pi * (
    R * R * R +
    (-2.) * K * K * R +
    K * (
      R * radical +
      K*K * math.atan2(R,radical)
    )
  )
