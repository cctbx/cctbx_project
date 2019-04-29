from __future__ import absolute_import, division, print_function
import math
pi = math.pi

def volume(radius_1, radius_2, distance):
  "http://mathworld.wolfram.com/Sphere-SphereIntersection.html"
  assert radius_1 >= 0
  assert radius_2 >= 0
  assert distance >= 0
  rrd = radius_1 + radius_2 - distance
  if (rrd <= 0):
    return 0
  if (distance <= abs(radius_1 - radius_2)):
    return 4/3 * pi * min(radius_1, radius_2)**3
  d = distance
  if (radius_1 == radius_2):
    return (pi * (4*radius_1 + d) * (2*radius_1 - d)**2) / 12
  if (radius_1 > radius_2):
    R, r = radius_1, radius_2
  else:
    R, r = radius_2, radius_1
  return pi * rrd**2 * (d**2 + 2*d*r - 3*r**2 + 2*d*R + 6*r*R - 3*R**2) \
       / (12*d)

def exercise():
  from libtbx.test_utils import approx_equal
  assert volume(0, 0, 0) == 0
  assert volume(1, 2, 4) == 0
  assert approx_equal(volume(4, 2, 1), 33.5103216383)
  assert approx_equal(volume(4, 4, 1), 218.078890037)
  assert approx_equal(volume(4+1.e-6, 4, 1), 218.078978001)
  assert approx_equal(volume(4, 4+1.e-6, 1), 218.078978001)
  assert approx_equal(volume(3, 4, 2), 94.9022780772)
  assert approx_equal(volume(4, 3, 2), 94.9022780772)
  assert approx_equal(volume(2, 2.3, 0), 33.5103216383)
  assert approx_equal(volume(2, 2.3, 0.3), 33.5103216383)
  assert approx_equal(volume(2, 2.3, 0.31), 33.505660942)
  assert approx_equal(volume(2, 2.3, 4.29), 0.000335811722418)
  assert approx_equal(volume(2, 2.3, 4.3), 0)
  print("OK")

if (__name__ == "__main__"):
  exercise()
