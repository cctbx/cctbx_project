from scitbx.math.ext import *
import scitbx.math.eigensystem
import scitbx.math.gaussian
from scitbx import matrix
from scitbx.array_family import flex
from stdlib import math

def euler_angles_as_matrix(angles, deg=00000):
  if (deg):
    angles = [a*math.pi/180 for a in angles]
  c = [math.cos(a) for a in angles]
  s = [math.sin(a) for a in angles]
  # Based on subroutine angro7 in CONVROT.
  # Urzhumtseva, L.M. & Urzhumtsev, A.G. (1997).
  return matrix.sqr((
     c[0]*c[1]*c[2]-s[0]*s[2],
    -c[0]*c[1]*s[2]-s[0]*c[2],
     c[0]*s[1],
     s[0]*c[1]*c[2]+c[0]*s[2],
    -s[0]*c[1]*s[2]+c[0]*c[2],
     s[0]*s[1],
    -s[1]*c[2],
     s[1]*s[2],
     c[1]))

class erf_verification:

  def __init__(self, tolerance=1.e-10):
    self.tolerance = tolerance
    self.max_delta = 0

  def __call__(self, f, x, expected_result):
    result = f(x)
    if (str(result) != str(expected_result)):
      me = str(result).lower().split("e")
      mex = str(expected_result).lower().split("e")
      assert len(me) == len(mex)
      if (len(me) == 2):
        assert me[1] == mex[1]
      delta = abs(float(me[0]) - float(mex[0]))
      if (self.max_delta < delta):
        self.max_delta = delta
        if (delta > self.tolerance):
          print x, expected_result, result, delta

class minimum_covering_sphere_nd:

  def __init__(self, points, epsilon):
    if (epsilon is None): epsilon = 1.e-6
    assert len(points) > 0
    n_dim = len(points[0].elems)
    w = 1./len(points)
    weights = matrix.row([w for i in xrange(len(points))])
    self._n_iterations = 0
    while 1:
      x = matrix.col([0]*n_dim)
      for w,t in zip(weights.elems,points):
        x += w * t
      radii = matrix.col([abs(x-t) for t in points])
      sigma = 0
      for w,r in zip(weights.elems,radii.elems):
        sigma += w * r**2
      sigma = math.sqrt(sigma)
      tau = radii.max()
      if (tau - sigma < tau * epsilon):
        break
      w_r = []
      for w,r in zip(weights.elems,radii.elems):
        w_r.append(w * r)
      w_r = matrix.col(w_r)
      sum_w_r = w_r.sum()
      assert sum_w_r != 0
      weights = w_r / sum_w_r
      self._n_iterations += 1
    self._center = x.elems
    self._radius = tau

  def n_iterations(self):
    return self._n_iterations

  def center(self):
    return self._center

  def radius(self):
    return self._radius

def minimum_covering_sphere(points, epsilon=None):
  if (epsilon is None): epsilon = 1.e-6
  if (isinstance(points, flex.vec3_double)):
    return minimum_covering_sphere_3d(points=points, epsilon=epsilon)
  else:
    return minimum_covering_sphere_nd(points=points, epsilon=epsilon)
