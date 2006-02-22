import scitbx.math
from scitbx import matrix
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal
from stdlib import math

def exercise(axis_range=2, angle_max_division=12, angle_min_power=-30):
  from_matrix = scitbx.math.r3_rotation_axis_and_angle_from_matrix(
    r=[1,0,0,0,1,0,0,0,1])
  assert approx_equal(from_matrix.axis, [1/3**0.5]*3)
  assert approx_equal(from_matrix.angle(deg=True), 0)
  assert approx_equal(from_matrix.angle(deg=False), 0)
  from_matrix = scitbx.math.r3_rotation_axis_and_angle_from_matrix(r=[0]*9)
  assert approx_equal(from_matrix.axis, [0,0,1])
  assert approx_equal(from_matrix.angle(deg=True), 90)
  assert approx_equal(from_matrix.angle(deg=False), math.pi/2)
  #
  angles = []
  for d in xrange(1,angle_max_division+1):
    angles.append(360/d)
    angles.append(-360/d)
  for p in xrange(-angle_min_power+1):
    angles.append(10**(-p))
    angles.append(-10**(-p))
  hex_orth = matrix.sqr([
    8.7903631196301042, -4.3951815598150503, 0,
    0, 7.6126777700894994, 0,
    0, 0, 14.943617303371177])
  for u in xrange(-axis_range, axis_range+1):
    for v in xrange(-axis_range, axis_range+1):
      for w in xrange(axis_range+1):
        for axis in [(u,v,w), (hex_orth*matrix.col((u,v,w))).elems]:
          for angle in angles:
            try:
              r = scitbx.math.r3_rotation_axis_and_angle_as_matrix(
                axis=axis, angle=angle, deg=True)
            except RuntimeError:
              assert axis == (0,0,0)
              break
            else:
              from_matrix = scitbx.math.r3_rotation_axis_and_angle_from_matrix(
                r=r)
              rr = from_matrix.as_matrix()
              assert approx_equal(rr, r)
              for deg in [False, True]:
                rr = scitbx.math.r3_rotation_axis_and_angle_as_matrix(
                  axis=from_matrix.axis,
                  angle=from_matrix.angle(deg=deg),
                  deg=deg,
                  min_axis_length=1-1.e-5)
                assert approx_equal(rr, r)
                try:
                  scitbx.math.r3_rotation_axis_and_angle_as_matrix(
                    axis=from_matrix.axis,
                    angle=from_matrix.angle(deg=deg),
                    deg=deg,
                    min_axis_length=1+1.e-5)
                except RuntimeError: pass
                else: raise RuntimeError("Exception expected.")
  #
  print format_cpu_times()

if (__name__ == "__main__"):
  exercise()
