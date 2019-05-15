from __future__ import absolute_import, division, print_function
from six.moves import range
def run(args):
  assert len(args) == 0
  n_trials = 100
  from scitbx.math.minimum_covering_ellipsoid import compute as mce_compute
  from scitbx.array_family import flex
  from libtbx.test_utils import approx_equal, is_below_limit
  # XXX point group 222 should be sufficient, but 432 is currently needed
  point_group_432_rotation_matrices = [
    (1,0,0,0,1,0,0,0,1),
    (1,0,0,0,0,-1,0,1,0),
    (1,0,0,0,0,1,0,-1,0),
    (0,0,1,0,1,0,-1,0,0),
    (0,0,-1,0,1,0,1,0,0),
    (0,-1,0,1,0,0,0,0,1),
    (0,1,0,-1,0,0,0,0,1),
    (0,0,1,1,0,0,0,1,0),
    (0,1,0,0,0,1,1,0,0),
    (0,-1,0,0,0,-1,1,0,0),
    (0,0,1,-1,0,0,0,-1,0),
    (0,-1,0,0,0,1,-1,0,0),
    (0,0,-1,-1,0,0,0,1,0),
    (0,0,-1,1,0,0,0,-1,0),
    (0,1,0,0,0,-1,-1,0,0),
    (1,0,0,0,-1,0,0,0,-1),
    (-1,0,0,0,1,0,0,0,-1),
    (-1,0,0,0,-1,0,0,0,1),
    (0,1,0,1,0,0,0,0,-1),
    (0,-1,0,-1,0,0,0,0,-1),
    (0,0,1,0,-1,0,1,0,0),
    (0,0,-1,0,-1,0,-1,0,0),
    (-1,0,0,0,0,1,0,1,0),
    (-1,0,0,0,0,-1,0,-1,0)]
  def check(center, radii, rotation):
    a,b,c = radii
    points_principal = flex.vec3_double([
      (-a,0,0),
      (a,0,0),
      (0,-b,0),
      (0,b,0),
      (0,0,-c),
      (0,0,c)])
    points = rotation * points_principal + center
    mce = mce_compute(points)
    assert approx_equal(mce.center, center)
    assert approx_equal(sorted(mce.radii), sorted(radii))
    assert approx_equal(mce.rotation.determinant(), 1)
    points_mce = mce.rotation.inverse().elems * (points - mce.center)
    rms = []
    for r in point_group_432_rotation_matrices:
      rp = r * points_mce
      rms.append(rp.rms_difference(points_principal))
    assert is_below_limit(value=min(rms), limit=1e-8, eps=0)
  mt = flex.mersenne_twister(seed=0)
  check((0,0,0), (1,2,3), (1,0,0,0,1,0,0,0,1))
  for i_trial in range(n_trials):
    center = list(mt.random_double(size=3)*8-4)
    radii = list(mt.random_double(size=3)*3+0.1)
    rotation = mt.random_double_r3_rotation_matrix()
    check(center, radii, rotation)
  from libtbx.utils import format_cpu_times
  print(format_cpu_times())

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
