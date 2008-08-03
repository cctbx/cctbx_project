from scitbx.array_family import flex
from scitbx.math import least_squares_plane
from scitbx import matrix as mat
from libtbx.utils import show_times_at_exit
import random
from libtbx.test_utils import approx_equal

def exercise():
  points = [ mat.col(x) for x in [(1, 2, 3), (-1, -2, -3), (1, 1, 1),
                                  (1.2, 2.1, 2.9), (-0.9, -2.1, -3.1),
                                  (1.1, 0.8, 1.2)] ]
  def distance(n,d):
    u = n/d
    return sum([ (u.dot(x) - 1)**2 for x in points ])
  flex_points = flex.vec3_double(points)
  p = least_squares_plane(flex_points, origin=(0,0,0))
  n = mat.col(p.normal)
  d = p.distance_to_origin
  assert approx_equal(abs(n), 1)
  dist0 = distance(n,d)
  for i in xrange(5000):
    d1 = d + random.uniform(-0.1, 0.1)
    n1 = mat.rec(flex.random_double_r3_rotation_matrix(), (3,3))*n
    dist = distance(n1, d1)
    assert dist > dist0

def run():
  show_times_at_exit()
  exercise()

if __name__ == '__main__':
  run()
