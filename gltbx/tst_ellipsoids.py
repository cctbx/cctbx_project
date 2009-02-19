import random
from scitbx.array_family import flex
import scitbx.math
from scitbx import matrix
from libtbx.test_utils import Exception_expected, approx_equal
from gltbx import quadrics

def exercise_ellipsoid(n_trials=100, n_sub_trials=10):
  rnd = random.Random(0)
  for i in xrange(n_trials):
    centre = matrix.col([ rnd.random() for k in xrange(3) ])
    half_lengths = matrix.col([ 0.1 + rnd.random() for k in xrange(3) ])
    r = scitbx.math.euler_angles_as_matrix(
      [ rnd.uniform(0, 360) for i in xrange(3) ], deg=True)
    metrics = r * matrix.diag([ x**2 for x in half_lengths ]) * r.transpose()
    t = quadrics.ellipsoid_to_sphere_transform(centre, metrics.as_sym_mat3())
    assert approx_equal(t.translation_part(), centre)
    m = matrix.sqr(t.linear_part())
    for j in xrange(n_sub_trials):
      y = matrix.col([ rnd.random() for k in xrange(3) ])
      c_y = y.transpose() * y
      x = m*y
      c_x = x.transpose() * metrics.inverse() * x
      assert approx_equal(c_x, c_y)
  r = scitbx.math.euler_angles_as_matrix((30, 115, 260), deg=True)
  centre = matrix.col((-1, 2, 3))
  metrics = r * matrix.diag((-1, 0.1, 1)) * r.transpose()
  t = quadrics.ellipsoid_to_sphere_transform(centre, metrics.as_sym_mat3())
  assert t.non_positive_definite()
  try:
    t.linear_part()
  except RuntimeError:
    pass
  else:
    raise Exception_expected

def run():
  exercise_ellipsoid()
  print "OK"


if __name__ == '__main__':
  run()
