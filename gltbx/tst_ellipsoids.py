#coding: utf-8
import random
from scitbx.array_family import flex
import scitbx.math
from scitbx import matrix
from libtbx.test_utils import approx_equal

def exercise_ellipsoid(n_trials=100, n_sub_trials=10):
  from gltbx import quadrics
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
    assert m.determinant() > 0
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
  x = r * matrix.col((1,0,0))
  assert x.transpose() * metrics.inverse() * x > 0

def time_ellipsoid(n=1000000):
  from gltbx.quadrics import time_ellipsoid_to_sphere_transform
  import timeit
  u = flex.sym_mat3_double(n, (0.0008, 0.0004, 0.0002,
                               0.0001, 0.00015, 0.00005))
  timer = timeit.Timer(lambda: time_ellipsoid_to_sphere_transform(u))
  print ("%i ellipsoid --> sphere transforms: %.3g μs per transform"
         % (n, timer.timeit(1)/n*1e6))

def run():
  import sys
  from libtbx.option_parser import option_parser
  try:
    import gltbx.gl
  except ImportError:
    print "Skipping gltbx/tst_ellipsoids.py: gltbx.gl module not available."
    sys.exit(1)

  exercise_ellipsoid()

  command_line = (option_parser(
    usage="",
    description="")
    .option(None, "--time", action="store_true")
    ).process(args=sys.argv[1:])
  if command_line.options.time:
    if command_line.args:
      time_ellipsoid(int(command_line.args[0]))
    else:
      time_ellipsoid()
  print "OK"

if __name__ == '__main__':
  run()
