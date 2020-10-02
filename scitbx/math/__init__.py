from __future__ import absolute_import, division, print_function

import functools
import math
import operator

from scitbx.math.ext import *
import scitbx.linalg.eigensystem
import scitbx.math.gaussian # implicit import
from scitbx import matrix
from scitbx.array_family import flex
from boost_adaptbx.boost import rational # implicit import
from six.moves import range
from six.moves import zip

gaussian_fit_1d_analytical.__doc__ = """
Fits a gaussian function to a list of points.

Parameters
----------
x : scitbx.array_family.flex.double
y : scitbx.array_family.flex.double

Returns
-------
scitbx.math.gaussian_fit_1d_analytical

Examples
--------
>>> from scitbx.array_family import flex
>>> from scitbx.math import gaussian_fit_1d_analytical
>>> x, y = flex.double(), flex.double()
>>> x.append(-1) ; y.append(.1)
>>> x.append(0) ; y.append(1)
>>> x.append(1) ; y.append(.1)
>>> fit = gaussian_fit_1d_analytical(x, y)
>>> print fit.a, fit.b
1.0 2.30258509299
"""

def median_statistics(data):
  return flex.median.dispersion(data)

class line_given_points(object):

  def __init__(self, points):
    self.points = [matrix.col(point) for point in points]
    self.delta = self.points[1] - self.points[0]
    self.delta_norm_sq = float(self.delta.norm_sq())

  def distance_sq(self, point):
    "http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html"
    if (self.delta_norm_sq == 0):
      return (point - self.points[0]).norm_sq()
    return self.delta.cross(point - self.points[0]).norm_sq() \
         / self.delta_norm_sq

  def perp_xyz(self, point):
    t = -(self.points[0] - point).dot(self.delta) / self.delta.norm_sq()
    return self.points[0]+(self.delta*t)

def r3_rotation_cos_rotation_angle_from_matrix(r):

  return max(
    min( ( r[0] + r[4] + r[8] - 1.0 ) / 2.0, 1.0 ),
    -1.0
    )

def r3_rotation_matrix_logarithm(r):
  """
  Chirikjian, G. Stochastic models, information theory, and Lie groups.
  Volume 2, pp 39-40.
  Applied and Numerical Harmonic Analysis, Birkhauser
  """

  cos_theta = r3_rotation_cos_rotation_angle_from_matrix( r )

  if cos_theta < 1.0:
      theta = math.acos( cos_theta )
      a1 = theta / math.sin( theta )

  else:
      a1 = 1.0

  return 0.5 * a1 * ( r - r.transpose() )

def r3_rotation_average_rotation_matrix_from_matrices(*matrices):
  """
  Calculates an approximation of the Riemannian (geometric) mean through
  quaternion averaging

  Arithmetic and geometric solutions for average rigid-body rotation (2010).
  I Sharf, A Wolf & M.B. Rubin, Mechanism and Machine Theory, 45, 1239-1251
  """

  if not matrices:
    raise TypeError("Average of empty sequence")

  import scitbx.matrix

  average = functools.reduce(
    operator.add,
    [ scitbx.matrix.col( r3_rotation_matrix_as_unit_quaternion( e ) )
      for e in matrices ]
    )

  return r3_rotation_unit_quaternion_as_matrix( average.normalize() )

def r3_rotation_average_rotation_via_lie_algebra(matrices, maxiter = 100, convergence = 1E-4):
  """
  Calculates an approximation of the Riemannian (geometric) mean by iteratively
  averaging the rotation matrix logarithms

  Nonparametric Second-order Theory of Error Propagation on Motion Groups (2008).
  Y. Wang & G. S. Chirikjian, The International Journal of Robotics Research 27, 1258-1273
  """

  if not matrices:
    raise TypeError("Average of empty sequence")

  if maxiter < 0 or convergence <= 0:
    raise ValueError("Invalid iteration parameters")

  from scitbx.math import so3_lie_algebra

  conv_sq = convergence * convergence
  norm = 1.0 / len( matrices )
  current = matrices[0]

  for i in range( maxiter ):
    inverse = current.inverse()
    log_diff = norm * functools.reduce(
      operator.add,
      [ so3_lie_algebra.element.from_rotation_matrix( matrix = inverse * m )
        for m in matrices ]
      )

    if log_diff.norm_sq() < conv_sq:
      return current

    current *= log_diff.exponential()

  raise RuntimeError("Iteration limit exceeded")

def euler_angles_as_matrix(angles, deg=False):
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

class erf_verification(object):

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
          print(x, expected_result, result, delta)

class minimum_covering_sphere_nd(object):

  def __init__(self,
        points,
        epsilon=None,
        radius_if_one_or_no_points=1,
        center_if_no_points=(0,0,0)):
    assert len(points) > 0 or radius_if_one_or_no_points >= 0
    if (epsilon is None): epsilon = 1.e-6
    self._n_iterations = 0
    if (len(points) == 0):
      self._center = center_if_no_points
      self._radius = radius_if_one_or_no_points
      return
    if (len(points) == 1):
      self._center = tuple(points[0])
      self._radius = radius_if_one_or_no_points
      return
    n_dim = len(points[0].elems)
    w = 1./len(points)
    weights = matrix.row([w for i in range(len(points))])
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

def minimum_covering_sphere(
      points,
      epsilon=None,
      radius_if_one_or_no_points=1,
      center_if_no_points=(0,0,0)):
  if (epsilon is None): epsilon = 1.e-6
  if (isinstance(points, flex.vec3_double)):
    impl = minimum_covering_sphere_3d
  else:
    impl = minimum_covering_sphere_nd
  return impl(
    points=points,
    epsilon=epsilon,
    radius_if_one_or_no_points=radius_if_one_or_no_points,
    center_if_no_points=center_if_no_points)

def row_echelon_back_substitution_int(
      row_echelon_form,
      v=None,
      solution=None,
      independent_flags=None):
  if (v is not None):
    assert v.nd() == 1
  return ext.row_echelon_back_substitution_int(
    row_echelon_form, v, solution, independent_flags)

def row_echelon_back_substitution_float(
      row_echelon_form,
      v=None,
      solution=None):
  assert solution is not None
  if (v is not None):
    assert v.nd() == 1
  return ext.row_echelon_back_substitution_float(
    row_echelon_form, v, solution)

def solve_a_x_eq_b_min_norm_given_a_sym_b_col(
      a, b,
      relative_min_abs_pivot=1e-12,
      absolute_min_abs_pivot=0,
      back_substitution_epsilon_factor=10):
  """\
Assumes a is symmetric, without checking to avoid overhead.

Special case of
  generalized_inverse(a) * b
taking advantage of the fact that a is real and symmetric.

Opportunistic algorithm: first assumes that a has full rank. If this
is true, solves a*x=b for x using simple back-substitution.

Only if a is rank-deficient:
  To obtain the x with minimum norm, transforms a to a basis formed
  by its eigenvectors, solves a*x=b in this basis, then transforms
  x back to the original basis system.

Returns None if a*x=b has no solution.
"""
  if (isinstance(a, matrix.rec)):
    a = a.as_flex_double_matrix()
  if (isinstance(b, matrix.rec)):
    assert b.n_columns() == 1
    b = flex.double(b)
  if (relative_min_abs_pivot is None):
    min_abs_pivot = 0
  else:
    min_abs_pivot = \
       max(a.all()) \
     * flex.max(flex.abs(a)) \
     * relative_min_abs_pivot
  min_abs_pivot = max(min_abs_pivot, absolute_min_abs_pivot)
  epsilon = min_abs_pivot * back_substitution_epsilon_factor
  aw = a.deep_copy()
  bw = b.deep_copy()
  ef = row_echelon_full_pivoting(
    a_work=aw, b_work=bw, min_abs_pivot=min_abs_pivot)
  if (ef.nullity == 0):
    x = ef.back_substitution(
      free_values=flex.double(ef.nullity), epsilon=epsilon)
  else:
    assert a.is_square_matrix()
    aw = a.deep_copy()
    es = scitbx.linalg.eigensystem.real_symmetric(aw)
    c = es.vectors() # may be left-handed, but that's OK here
    ct = c.matrix_transpose()
    aw = c.matrix_multiply(aw).matrix_multiply(ct)
    bw = c.matrix_multiply(b)
    max_rank = ef.rank
    ef = row_echelon_full_pivoting(a_work=aw, b_work=bw, max_rank=max_rank)
    assert ef.rank == max_rank
    x = ef.back_substitution(
      free_values=flex.double(ef.nullity, 0), epsilon=epsilon)
    if (x is not None):
      x = ct.matrix_multiply(x)
  return x

class finite_difference_computation(object):
  """ Computing the derivative of a scalar function f(x, ...) wrt x,
  where x may be a vector, involves computing f(x + h, ...)
  where h is suitably chosen a delta. So as to keep numerical accuracy, it
  is best to follow a well known recipe to choose h knowing its
  direction. When using that value, a finite difference computation typically
  reaches a precision of self.precision

  Reference: chap 5.7 of Numerical Recipes in C
  """

  def __init__(self):
    e = math.sqrt(scitbx.math.double_numeric_limits.epsilon)
    self.precision = math.sqrt(e)

  def best_delta(self, point, direction):
    """ Compute the value h in the given direction that is most suitable
    for use in finite difference computation around the given point x.
    The vector 'direction' may not be a unit-vector.
    """
    x, s = point, self.precision
    u = direction/abs(direction)
    h = s*u
    h = (x + h) - x
    return h

def point_on_sphere(r, s_deg, t_deg, center):
  """
  Returns a point on sphere of radius r: see definition of polar (sperical)
  coordinate system in 3D.
  """
  sin = math.sin
  cos = math.cos
  s = s_deg*math.pi/180
  t = t_deg*math.pi/180
  x = center[0] + r * cos(s) * sin(t)
  y = center[1] + r * sin(s) * sin(t)
  z = center[2] + r * cos(t)
  return x,y,z

def equally_spaced_points_on_vector(start, end, n=None, step=None):
  assert [n, step].count(None) == 1
  vec = [end[0]-start[0],end[1]-start[1],end[2]-start[2]]
  r = flex.vec3_double([vec])
  if(n is not None):
    assert n > 0
  else:
    assert step > 0
    vec_length = math.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)
    n = int(vec_length/step)-1
  dr = r*(1/float(n+1))
  points = flex.vec3_double()
  for i in range(n+1):
    points.extend(dr * i + start)
  points.append(end)
  return points

def five_number_summary(data):
  """
  Returns the Tukey five number summary (min, lower hinge, median, upper hinge,
  max) for a sequence of observations. This function gives the same results
  as R's fivenum function.
  """
  try:
    sorts = flex.sorted(data)
  except AttributeError:
    sorts = sorted(data)
  n = len(sorts)
  if n % 2:
    med = sorts[n // 2]
    lower = sorts[:((n // 2) + 1)]
    upper = sorts[(n // 2):]
  else:
    med = (sorts[n//2] + sorts[n//2 - 1]) / 2
    lower = sorts[:(n // 2)]
    upper = sorts[(n // 2):]
  n = len(lower)
  if n % 2:
    lhinge = lower[n // 2]
    uhinge = upper[n // 2]
  else:
    lhinge = (lower[n//2] + lower[n//2 - 1]) / 2
    uhinge = (upper[n//2] + upper[n//2 - 1]) / 2
  return sorts[0], lhinge, med, uhinge, sorts[-1]

class sin_cos_table(object):
  """
  Always evaludate performance when it is important!
  """
  def __init__(self, n):
    self.n = n
    self.step = 2*math.pi/self.n
    self.sin_table = flex.double()
    self.cos_table = flex.double()
    for i in range(self.n):
      self.sin_table.append(math.sin(i*self.step))
      self.cos_table.append(math.cos(i*self.step))

def similarity_indices(x, eps):
  """
  Best explained by examples:
    for array [0.232, 0.002, 0.45, 1.2, 0.233, 1.2, 0.5, 0.231, 0,0,0,0,0]
    return    [0,     1,     2,    3,   0,     3,   4,   0,     1,1,1,1,1]
    if eps=0.01.
  Other examples:
    [0, 2, 0, 2]
    [0, 1, 0, 1]

    [0, 0, 0, 0]
    [0, 0, 0, 0]

    [1, 2, 3, 4]
    [0, 1, 2, 3]
  """
  ii_all = []
  for i, xi in enumerate(x):
    ii = []
    for j, xj in enumerate(x):
      if(abs(xi-xj)<eps): ii.append(j)
    if(not ii in ii_all): ii_all.append(ii)
  y = x[:]
  for i, ii in enumerate(ii_all):
    for iii in ii: y[iii]=i
  return y

def simpson(f, a, b, n):
  """
  Integral of f on [a,b] using composite Simpson's rule.
  """
  if(n % 2 != 0): raise ValueError("n must be even" % n)
  h = (b - a) / n
  s = f(a) + f(b)
  for i in range(1, n, 1):
    if(i%2==0): s += 4 * f(a + i * h)
    else:       s += 2 * f(a + i * h)
  return s * h / 3
