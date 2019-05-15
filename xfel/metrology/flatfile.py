# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

# TODO (1) OUTPUT TO SENSIBLE COORDINATES (2) WEIGHTING


"""
For Mikhail:

  Quad 3, pair 0 should have S2 20907 (not 20904).

  Quad 3, pair 4 should have L2 43477 (not 43476), but that may be
  rounding error.

  Quad 3, pair 0 has L2 43857 which is 300 um longer than others

  XXX Really need to be able to visualise this in 3D!

Other open issues:

  Possible to participate in next optical measurement?

  Is the center of the rectangle the center of the pixel arrays?
"""

from __future__ import absolute_import, division, print_function

import copy
import math
import sys

from scitbx import matrix
from scitbx.array_family import flex

#from scitbx import matrix
#from scitbx.array_family import flex
#from xfel.cxi.cspad_ana.parse_calib import calib2sections

from scitbx import lbfgs
from six.moves import range


class optimise_rectangle(object):
  """Misnomer, optimises everything."""

  def __init__(self, vertices, weights, c, u, v, l):
    """Not really, it's actually first and second coordinate in
    whatever coordinate frame is in use, most likely x and y.  First
    and second spatial coordinate.  XXX Try lbfgsb, for constrained
    optimisation?

    @param coords Data to optimise, five-tuple
    @param theta  Starting guess, rotation angle
    @param t_s    Starting guess, slow displacement
    @param t_f    Starting guess, fast displacement
    """

    self.x = flex.double((c(0, 0), c(1, 0), c(2, 0),
                          u(0, 0), u(1, 0), u(2, 0),
                          v(0, 0), v(1, 0), v(2, 0)))
    self.vertices = vertices
    self.weights = weights
    self.minimizer = lbfgs.run(target_evaluator=self)


  def compute_functional_and_gradients(self):
    """Compute sum of squared residual, and derivatives of theta,
    slow, and fast.  XXX THIS IS ALL DIFFERENT FROM MY PEN-AND-PAPER
    WORK"""

    c = matrix.col((self.x[0], self.x[1], self.x[2]))
    u = matrix.col((self.x[3], self.x[4], self.x[5]))
    v = matrix.col((self.x[6], self.x[7], self.x[8]))

    f =   self.weights[0] * (self.vertices[0] - c + v).norm_sq() \
        + self.weights[1] * (self.vertices[1] - c - u).norm_sq() \
        + self.weights[2] * (self.vertices[2] - c - v).norm_sq() \
        + self.weights[3] * (self.vertices[3] - c + u).norm_sq() \
        + (u.norm_sq() - v.norm_sq())**2

    dfdc = -2 * self.weights[0] * (self.vertices[0] - c + v) \
        -   2 * self.weights[1] * (self.vertices[1] - c - u) \
        -   2 * self.weights[2] * (self.vertices[2] - c - v) \
        -   2 * self.weights[3] * (self.vertices[3] - c + u)
    dfdu = -2 * self.weights[1] * (self.vertices[1] - c - u) \
        +   2 * self.weights[3] * (self.vertices[3] - c + u) \
        +   4 * (u.norm_sq() - v.norm_sq()) * u
    dfdv = +2 * self.weights[0] * (self.vertices[0] - c + v) \
        -   2 * self.weights[2] * (self.vertices[2] - c - v) \
        -   4 * (u.norm_sq() - v.norm_sq()) * v

    g = flex.double((dfdc(0, 0), dfdc(1, 0), dfdc(2, 0),
                     dfdu(0, 0), dfdu(1, 0), dfdu(2, 0),
                     dfdv(0, 0), dfdv(1, 0), dfdv(2, 0)))
    return (f, g)


class optimise_rectangle_2(object):
  """Misnomer, optimises everything."""

  def __init__(self, vertices, weights, c, u, v, l):
    """Not really, it's actually first and second coordinate in
    whatever coordinate frame is in use, most likely x and y.  First
    and second spatial coordinate.  XXX Try lbfgsb, for constrained
    optimisation?

    @param coords Data to optimise, five-tuple
    @param theta  Starting guess, rotation angle
    @param t_s    Starting guess, slow displacement
    @param t_f    Starting guess, fast displacement
    """

    self.x = flex.double((c(0, 0), c(1, 0), c(2, 0),
                          u(0, 0), u(1, 0), u(2, 0),
                          v(0, 0), v(1, 0), v(2, 0)))
    self.vertices = vertices
    self.weights = weights
    self.minimizer = lbfgs.run(target_evaluator=self)


  def compute_functional_and_gradients(self):
    """Compute sum of squared residual, and derivatives of theta,
    slow, and fast.  XXX THIS IS ALL DIFFERENT FROM MY PEN-AND-PAPER
    WORK"""

    c = matrix.col((self.x[0], self.x[1], self.x[2]))
    u = matrix.col((self.x[3], self.x[4], self.x[5]))
    v = matrix.col((self.x[6], self.x[7], self.x[8]))

    f =   self.weights[0] * (self.vertices[0] - c + v).length() \
        + self.weights[1] * (self.vertices[1] - c - u).length() \
        + self.weights[2] * (self.vertices[2] - c - v).length() \
        + self.weights[3] * (self.vertices[3] - c + u).length() \
        + (u.norm_sq() - v.norm_sq())**2

    # XXX normalize() will fail for zero-vectors.
    dfdc = -1 * self.weights[0] * (self.vertices[0] - c + v).normalize() \
        -   1 * self.weights[1] * (self.vertices[1] - c - u).normalize() \
        -   1 * self.weights[2] * (self.vertices[2] - c - v).normalize() \
        -   1 * self.weights[3] * (self.vertices[3] - c + u).normalize()
    dfdu = -1 * self.weights[1] * (self.vertices[1] - c - u).normalize() \
        +   1 * self.weights[3] * (self.vertices[3] - c + u).normalize() \
        +   4 * (u.norm_sq() - v.norm_sq()) * u
    dfdv = +1 * self.weights[0] * (self.vertices[0] - c + v).normalize() \
        -   1 * self.weights[2] * (self.vertices[2] - c - v).normalize() \
        -   4 * (u.norm_sq() - v.norm_sq()) * v

    g = flex.double((dfdc(0, 0), dfdc(1, 0), dfdc(2, 0),
                     dfdu(0, 0), dfdu(1, 0), dfdu(2, 0),
                     dfdv(0, 0), dfdv(1, 0), dfdv(2, 0)))
    return (f, g)


class optimise_rectangle_3(object):

  def __init__(self, vertices, weights, u, v, side_long, side_short):
    """
    Four 3D vertices, four weights.  u and v are 3D vectors.
    Optimises weights, u, and v.

    @param vertices Must be ordered in counter-clockwise direction
    @param weights  Four scalars corresponding to @p vertices XXX or w?
    @param u        3D vector
    @param v        3D vector
    """

    self.x = flex.double(
      (weights.elems[0], weights.elems[1], weights.elems[2], weights.elems[3],
       u.elems[0], u.elems[1], u.elems[2],
       v.elems[0], v.elems[1], v.elems[2]))

    self._vertices = vertices
    self._weights = weights

    self._side_long = side_long
    self._side_short = side_short

    self.minimizer = lbfgs.run(target_evaluator=self)


  def compute_functional_and_gradients(self):

    """
    while True:
      t = self.x[0]**2 + self.x[1]**2 + self.x[2]**2 + self.x[3]**2
#      print "Got weights", self.x[0]**2 / t, self.x[1]**2 / t, self.x[2]**2 / t, self.x[3]**2 / t

      for i in range(4):
        w = self.x[i]**2 / t

        ok = True
        self.x[i] = math.sqrt(w)
        if w > 0.35:
          self.x[i] = 0.3 / w * math.sqrt(w)
          ok = False

        if ok:
          break
    """

    w1 = flex.double((self.x[0],
                      self.x[1],
                      self.x[2],
                      self.x[3]))

    w2 = flex.double((w1[0]**2,
                      w1[1]**2,
                      w1[2]**2,
                      w1[3]**2))

    w2_tot = w2[0] + w2[1] + w2[2] + w2[3]
    assert w2_tot > 0 # XXX Pathology
    # w2 /= w2_tot # XXX Dangerous!  Don't do this!

    u = matrix.col((self.x[4], self.x[5], self.x[6]))
    v = matrix.col((self.x[7], self.x[8], self.x[9]))

    # Compute optimal rectangle center, with weights.  NOTE: center is
    # a function of (p, u, v, w).  See pen-and-paper derivation,
    # obtained from df/dc = 0.
    c = (w2[0] * (self._vertices[0] + v) +
         w2[1] * (self._vertices[1] - u) +
         w2[2] * (self._vertices[2] - v) +
         w2[3] * (self._vertices[3] + u)) / w2_tot

    # Weighted residual, the function (XXX functional) to optimize.
    f = (w2[0] * (self._vertices[0] - c + v).norm_sq() +
         w2[1] * (self._vertices[1] - c - u).norm_sq() +
         w2[2] * (self._vertices[2] - c - v).norm_sq() +
         w2[3] * (self._vertices[3] - c + u).norm_sq()) / w2_tot
    #f += (u.norm_sq() - v.norm_sq())**2
    #f += math.fabs(u.norm_sq() - v.norm_sq())
    f += (u.length() - v.length())**2
    f += ((u + v).length() - self._side_long)**2
    f += ((u - v).length() - self._side_short)**2
#    f += math.fabs((u + v).norm_sq() - self._side_long**2)
#    f += math.fabs((u - v).norm_sq() - self._side_short**2)

#    f += math.fabs(u.length() - v.length())
#    f += math.fabs((u + v).length() - self._side_long)
#    f += math.fabs((u - v).length() - self._side_short)

#    print "*** CHECK", math.fabs((u - v).length() - self._side_short), (u - v).length(), self._side_short

    # Derivatives of f w.r.t. weights
    t = 2 / w2_tot * (w2[0] * (self._vertices[0] - c + v) +
                      w2[1] * (self._vertices[1] - c - u) +
                      w2[2] * (self._vertices[2] - c - v) +
                      w2[3] * (self._vertices[3] - c + u))

    dfdw0 = 2 * w1[0] / w2_tot * (
      (self._vertices[0] - c + v).norm_sq() - t.dot(self._vertices[0] + v) + f)
    dfdw1 = 2 * w1[1] / w2_tot * (
      (self._vertices[1] - c - u).norm_sq() - t.dot(self._vertices[1] - u) + f)
    dfdw2 = 2 * w1[2] / w2_tot * (
      (self._vertices[2] - c - v).norm_sq() - t.dot(self._vertices[2] - v) + f)
    dfdw3 = 2 * w1[3] / w2_tot * (
      (self._vertices[3] - c + u).norm_sq() - t.dot(self._vertices[3] + u) + f)

    # These are vectors.  XXX Need to divide these by w2_tot
    dfdu = -2 * (w2[0] * (self._vertices[0] - c + v) +
                 w2[2] * (self._vertices[2] - c - v)) * (w2[3] - w2[1]) / w2_tot \
                 - 2 * w2[1] * (self._vertices[1] - c - u) * (w2[0] + w2[2] + 2 * w2[3]) / w2_tot \
                 + 2 * w2[3] * (self._vertices[3] - c + u) * (w2[0] + 2 * w2[1] + w2[2]) / w2_tot
    dfdu = dfdu / w2_tot
    dfdu += 2 * (1 - v.length() / u.length()) * u
    dfdu += 2 * (1 - self._side_long / (u + v).length()) * (u + v)
    dfdu += 2 * (1 - self._side_short / (u - v).length()) * (u - v)
#    dfdu += 4 * (u.norm_sq() - v.norm_sq()) * u
#    dfdu += +math.copysign(2, u.norm_sq() - v.norm_sq()) * u
#    dfdu += +math.copysign(2, (u + v).norm_sq() - self._side_long**2) * (u + v)
#    dfdu += +math.copysign(2, (u - v).norm_sq() - self._side_short**2) * (u - v)

#    if math.fabs(u.length() - v.length()) > 1e-6:
#      dfdu += (1 - v.length() / u.length()) / math.fabs(u.length() - v.length()) * u
#    dfdu += (1 - self._side_long / (u + v).length()) / math.fabs((u + v).length() - self._side_long) * (u + v)
#    dfdu += (1 - self._side_short / (u - v).length()) / math.fabs((u - v).length() - self._side_short) * (u - v)


    dfdv = -2 * (w2[1] * (self._vertices[1] - c - u) +
                 w2[3] * (self._vertices[3] - c + u)) * (w2[0] - w2[2]) / w2_tot \
                 + 2 * w2[0] * (self._vertices[0] - c + v) * (w2[1] + 2 * w2[2] + w2[3]) / w2_tot \
                 - 2 * w2[2] * (self._vertices[2] - c - v) * (2 * w2[0] + w2[1] + w2[3]) / w2_tot
    dfdv = dfdv / w2_tot
    dfdv +=  2 * (1 - u.length() / v.length()) * v
    dfdv +=  2 * (1 - self._side_long / (u + v).length()) * (u + v)
    dfdv += -2 * (1 - self._side_short / (u - v).length()) * (u - v)
#    dfdv += -4 * (u.norm_sq() - v.norm_sq()) * v
#    dfdv += -math.copysign(2, u.norm_sq() - v.norm_sq()) * v
#    dfdv += +math.copysign(2, (u + v).norm_sq() - self._side_long**2) * (u + v)
#    dfdv += -math.copysign(2, (u - v).norm_sq() - self._side_short**2) * (u - v)

#    if math.fabs(u.length() - v.length()) > 1e-6:
#      dfdv += (1 - u.length() / v.length()) / math.fabs(u.length() - v.length()) * v
#    dfdv += (1 - self._side_long / (u + v).length()) / math.fabs((u + v).length() - self._side_long) * (u + v)
#    dfdv += (1 - self._side_short / (u - v).length()) / math.fabs((u - v).length() - self._side_short) * (u - v) * (-1)

    g = flex.double(
      (dfdw0, dfdw1, dfdw2, dfdw3,
       dfdu.elems[0], dfdu.elems[1], dfdu.elems[2],
       dfdv.elems[0], dfdv.elems[1], dfdv.elems[2]))

    #print "*** DID ITERATION ***"

    return (f, g)


def _find_long_short(quadrants,
                     l_int=(-float('inf'), +float('inf')),
                     s_int=(-float('inf'), +float('inf'))):
  """The _find_long_short() function returns the average length of the
  long and short sides of the rectangles whose vertices are given in
  @p quadrants.

  @param quadrants XXX
  @param l_int     XXX Permitted interval for long (XXX even?) side
  @param s_int     XXX Permitted interval for short (XXX odd?) side
  @return          Three-tuple of the long side length, short side
                   length, and standard deviation from the mean XXX Could really return a struct (erm... class instance)
  """

  # Now assumes correct ordering.  XXX Does this do adequate checking?
  l_sum = l_ssq = 0
  s_sum = s_ssq = 0
  l_N = 0
  s_N = 0

  for (q, sensors) in quadrants.iteritems():
    for (s, vertices) in sensors.iteritems():
      if len(vertices) != 4:
        raise RuntimeError("Sensor in quadrant does not have four vertices")
      l = (vertices[1] - vertices[0]).length()
#      print "SIDE long:", l
      if l > l_int[0] and l < l_int[1]:
        l_sum += l
        l_ssq += l**2
        l_N += 1

      l = (vertices[3] - vertices[2]).length()
#      print "SIDE long:", l
      if l > l_int[0] and l < l_int[1]:
        l_sum += l
        l_ssq += l**2
        l_N += 1

      s = (vertices[2] - vertices[1]).length()
#      print "SIDE short:", s
      if s > s_int[0] and s < s_int[1]:
        s_sum += s
        s_ssq += s**2
        s_N += 1

      s = (vertices[0] - vertices[3]).length()
#      print "SIDE short:", s
      if s > s_int[0] and s < s_int[1]:
        s_sum += s
        s_ssq += s**2
        s_N += 1

  l_var = (l_ssq - l_sum**2 / l_N) / (l_N - 1)
  s_var = (s_ssq - s_sum**2 / s_N) / (s_N - 1)
  l_mu = l_sum / l_N
  s_mu = s_sum / s_N
  return (l_mu, math.sqrt(l_var), l_N, s_mu, math.sqrt(s_var), s_N)



  # Compute all the side lengths from all the sensors in all quadrants
  # using every pair of successive vertices.
  d = []
  for (q, sensors) in quadrants.iteritems():
    for (s, vertices) in sensors.iteritems():
      if len(vertices) != 4:
        raise RuntimeError("Sensor in quadrant does not have four vertices")
      for i in range(len(vertices)):
        d.append((vertices[i] - vertices[(i + 1) % len(vertices)]).length())
  if len(d) < 4:
    return (0, 0, 0)

  # Sort the side-lengths in-place, and compute separate averages for
  # the lower and upper half of the array.
  d.sort()
  h = len(d) // 2

  l_mu = s_mu = 0
  for i in range(h):
    s_mu += d[i]
  s_mu /= h
  for i in range(h, 2 * h):
    l_mu += d[i]
  l_mu /= h

  # Compute standard deviation XXX biased or unbiased?
  sigma = 0
  for i in range(h):
    sigma += (d[i] - s_mu)**2
  for i in range(h, 2 * h):
    sigma += (d[i] - l_mu)**2
  sigma = math.sqrt(sigma / (2 * h - 2))

  return (l_mu, s_mu, sigma)


def _find_aspect_ratio(quadrants, r_int=(-float('inf'), +float('inf'))):
  """The _find_aspect_ratio() function returns the average aspect
  ratio (long side divided by short side) of the rectangles whose
  vertices are given in @p quadrants.

  @param quadrants XXX
  @param r_int     XXX Permitted interval for aspect ratio
  @return          Three-tuple of the long side length, short side
                   length, and standard deviation from the mean XXX Could really return a struct (erm... class instance)
  """

  # Now assumes correct ordering.  XXX Does this do adequate checking?
  r_sum = r_ssq = 0
  r_N = 0

  for (q, sensors) in quadrants.iteritems():
    for (s, vertices) in sensors.iteritems():
      if len(vertices) != 4:
        raise RuntimeError("Sensor in quadrant does not have four vertices")
      s = (vertices[2] - vertices[1]).length() + (vertices[0] - vertices[3]).length()
      if s > 0:
        l = (vertices[1] - vertices[0]).length() + (vertices[3] - vertices[2]).length()
        r = l / s
        if r > r_int[0] and r < r_int[1]:
          r_sum += r
          r_ssq += r**2
          r_N += 1

  r_var = (r_ssq - r_sum**2 / r_N) / (r_N - 1)
  r_mu = r_sum / r_N
  return (r_mu, math.sqrt(r_var), r_N)


def _order_sensors(sensors):
  """The _order_senors() function reshuffles the sensors and vertices
  within a quadrant to match the indexing in the XTC stream.  The
  coordinates themselves are not changed.

  The _longest_side() function returns the coordinate index of the
  longest side of the polygon defined by the vertices in @p.  Here, it
  should be either 0 for x or 1 for y.

  The _check_rectangle() function verifies that the vertices, when
  traversed in order of increasing indices define normals that point
  along the positive z-axis.  XXX Name?!

  Returns dictionary

  @param sensors List of 4-tuples of 3D-vertices of a polygon
  @return  @c True if the normals of all vertices point along the
           positive z-axis

  """

  from collections import deque

  # Split the set of sensors into two disjoint (XXX) subsets depending
  # on orientation: either x-major or y-major.  Record their centers
  # for subsequent sorting.
  sensors_x = []
  sensors_y = []
  for v in sensors:
    if len(v) == 0:
      continue

    # Compute the center of the sensor.
    m = v[0]
    for j in range(1, len(v)):
      m += v[j]
    m /= len(v)

    # Order the vertices in counter-clockwise direction in the
    # xy-plane.
    sorted(v, key=lambda x: math.atan2((x - m).elems[1], (x - m).elems[0]))

    # Compute extent in all three dimensions, and determine longest
    # dimension.  l = 0, sensor is oriented along the x-axis, l = 1,
    # sensor is oriented along the y-axis.
    l = flex.double()
    for i in range(len(m.elems)):
      c = [v[j].elems[i] for j in range(len(v))]
      l.append(max(c) - min(c))
    l = flex.max_index(l)

    if l == 0:
      sensors_x.append((v, m))
    elif l == 1:
      sensors_y.append((v, m))
    else:
      raise RuntimeError("XXX This should not happen")

  # Sort x-major sensors by increasing y-value of the center, and
  # y-major sensors by increasing x-value of the center.
  sorted(sensors_x, key=lambda v_m: v_m[1].elems[1])
  sorted(sensors_y, key=lambda v_m1: v_m1[1].elems[0])

  # Write a dictionary with the sensors in the correct order, and
  # discard the sensor centers and convert to deque.  Now remains only
  # to cyclically permute the vertices within the sensors.
  #
  # 5 4 6 6  x
  # 5 4 7 7  ^
  # 2 2 0 1  |
  # 3 3 0 1  |
  #     y <--+
  assert len(sensors_x) == 4 and len(sensors_y) == 4 # XXX
  sensors_dict = {
    0: deque(sensors_x[1][0]), 1: deque(sensors_x[0][0]),
    2: deque(sensors_y[1][0]), 3: deque(sensors_y[0][0]),
    4: deque(sensors_x[2][0]), 5: deque(sensors_x[3][0]),
    6: deque(sensors_y[3][0]), 7: deque(sensors_y[2][0])}

  # Reshuffle vertices in sensors by cyclic permutation.
  for (i, v) in sensors_dict.iteritems():
    assert len(v) >= 3 # XXX
    d1 = v[1] - v[0]
    d2 = v[2] - v[1]

    # Ensure v[0] starts a long edge.
    if d1.length() < d2.length():
      d1 = d2
      v.rotate(-1)

    # Ensure v[0] corresponds to origin.  XXX Explain better.
    if    ((i == 0 or i == 1) and d1.elems[0] < 0) or \
          ((i == 2 or i == 3) and d1.elems[1] > 0) or \
          ((i == 4 or i == 5) and d1.elems[0] > 0) or \
          ((i == 6 or i == 7) and d1.elems[1] < 0):
      v.rotate(2)
    sensors_dict[i] = list(v)

  return sensors_dict


def _find_weight(p0, p1, p2, l, s):
  #print "Finding weight using %f and %f" % (l, s)


  # Rotate the p0 and p2 into the z-plane, use p1 as origin Find
  # rotation axis: this is the cross product of n and [0, 0, 1], where
  # n is the cross product of (p0 - p1) and (p2 - p1).  Rotation angle
  # is the arccos of the dot product of n and [0, 0, 1].
  n = (p0 - p1).cross(other=(p2 - p1)).normalize()
  r = matrix.col((+n(1, 0), -n(0, 0), 0))
  if r.length() > 0:
    R = r.axis_and_angle_as_r3_rotation_matrix(angle=math.acos(n(2, 0)))
  else:
    print("TOOK funny branch")
    R = matrix.identity(3)

  # Rotate, and make p1 the origin.  This is now a 2D-problem.
  q0 = matrix.col((R(0, 0) * (p0 - p1)(0, 0) + R(0, 1) * (p0 - p1)(1, 0),
                   R(1, 0) * (p0 - p1)(0, 0) + R(1, 1) * (p0 - p1)(1, 0)))
  q1 = matrix.col((0, 0))
  q2 = matrix.col((R(0, 0) * (p2 - p1)(0, 0) + R(0, 1) * (p2 - p1)(1, 0),
                   R(1, 0) * (p2 - p1)(0, 0) + R(1, 1) * (p2 - p1)(1, 0)))

  #if abs(q0(2, 0)) > 1e-7 or abs(q2(2, 0)) > 1e-7:
  #  print "B0RKED", q0(2, 0), q2(2, 0)


  # Find best fit against [l, 0] and [0, s]
  theta = math.atan2(
    q0(1, 0) * l - q2(0, 0) * s, q0(0, 0) * l + q2(1, 0) * s)
  r1 =  (q0 - l * matrix.col((+math.cos(theta), +math.sin(theta)))).length() \
      + (q2 - s * matrix.col((-math.sin(theta), +math.cos(theta)))).length()

  theta = math.atan2(
    q0(1, 0) * s - q2(0, 0) * l, q0(0, 0) * s + q2(1, 0) * l)
  r2 =  (q0 - s * matrix.col((+math.cos(theta), +math.sin(theta)))).length() \
      + (q2 - l * matrix.col((-math.sin(theta), +math.cos(theta)))).length()


  r_min = min(r1, r2)
  r_max = max(r1, r2)

  #print "Got", math.sqrt(r_min), math.sqrt(r_max)

  # Alternative: find best fit against [s, 0] and [0, l]

  # Return best fit of the two

  return math.sqrt(r_min)
#  return 1
#  return 1 / (1 + math.sqrt(r_min))


def _fit_rectangle(vertices, side_long, side_short): #p1, p2, p3, p4):
  """Least-squares fit a rectangle to four points.  Points must be
  given in counter-clockwise order.  XXX Check again all against
  pen-and-paper!

  XXX Need three vertices for plane.

  XXX IMPORTANT: the rectangle size is not constrained?!
  """

  from sys import float_info

  """

  #Mikhail's statistics, from
  # https://confluence.slac.stanford.edu/display/PCDS/2011-08-10+CSPad+alignment+parameters
  # Only minor discrepancies (typos?), and reordering of S1/S2, L1/L2,
  # etc, and signs.  Assumes the z-coordinate has been discarded XXX
  # do it here instead?  Make this it's own function?

  # Discard the z-component (XXX coordinate?)
  v = []
  for vertex in vertices:
    v.append(matrix.col((vertex(0, 0), vertex(1, 0))))

  L1 = (v[1] - v[0]).length()
  L2 = (v[3] - v[2]).length()
  S1 = (v[2] - v[1]).length()
  S2 = (v[0] - v[3]).length()

  if abs((v[1] - v[0])(0, 0)) > abs((v[1] - v[0])(1, 0)):
    # For horizontal sensors
    orientation = "H"
    dL1 = v[1][0] - v[2][0]
    dL2 = v[0][0] - v[3][0]
    dS2 = v[1][1] - v[0][1]
    dS1 = v[2][1] - v[3][1]
  else:
    # For vertical sensors
    orientation = "V"
    dL1 = v[1][1] - v[2][1]
    dL2 = v[0][1] - v[3][1]
    dS1 = v[1][0] - v[0][0]
    dS2 = v[2][0] - v[3][0]

  dS_L = (dS1 + dS2) / (L1 + L2)
  angle = math.atan(dS_L) * 180 / math.pi
  print "%c: S1 % 6d S2 % 6d dS1 % 6d dS2 % 6d " \
      "L1 % 6d L2 % 6d dL1 % 6d dL2 % 6d " \
      "<dS/L> % 8.5f angle % 8.5f"% \
      (orientation,
       int(S1), int(S2), int(dS1), int(dS2),
       int(L1), int(L2), int(dL1), int(dL2),
       dS_L, angle)
  return
  """

  # The center is always the mean of the four corners.
  c = 0.25 * (vertices[0] + vertices[1] + vertices[2] + vertices[3])

  # Determine Lagrange multiplier, lambda.
  d02 = (vertices[0] - vertices[2]).norm_sq()
  d13 = (vertices[1] - vertices[3]).norm_sq()

  # XXX Establish absolute tolerance for equality.  Square root of
  # smallest increment from maximum.  XXX Use something from cctbx
  # instead?
  tol = math.sqrt(math.ldexp(
      float_info.epsilon,
      math.frexp(flex.max(flex.abs(flex.double([d02, d13]))))[1] - 1))

  if abs(d02 - d13) < tol:
    l = 0
  else:
    l1 = 2 / (d02 - d13) * (-(d02 + d13) + 2 * math.sqrt(d02 * d13))
    l2 = 2 / (d02 - d13) * (-(d02 + d13) - 2 * math.sqrt(d02 * d13))
    if abs(l1) < 2 and abs(l2) >= 2:
      l = l1
    elif abs(l1) >= 2 and abs(l2) < 2:
      l = l2
    else:
      raise RuntimeError("Cannot determine Lagrange multiplier")

  u = (vertices[1] - vertices[3]) / (2 + l)
  v = (vertices[2] - vertices[0]) / (2 - l)
  r =   (vertices[0] - c + v).norm_sq() \
      + (vertices[1] - c - u).norm_sq() \
      + (vertices[2] - c - v).norm_sq() \
      + (vertices[3] - c + u).norm_sq() \
      + l * (u.norm_sq() - v.norm_sq())

  w = (u - v).length()
  h = (u + v).length()
  if w < h:
    t = h
    h = w
    w = t

#  print "%7.3f %7.3f %7.3f %7.3f -- %d %d" % (
#    (vertices[0] - c + v).length(),
#    (vertices[1] - c - u).length(),
#    (vertices[2] - c - v).length(),
#    (vertices[3] - c + u).length(),
#    w, h)


  # WEIGHT DETERMINATION, use (side_long, side_short) or (w, h) --
  # USELESS?

  weights = []
  for i in range(len(vertices)):
    weights.append(_find_weight(vertices[(i - 1) % len(vertices)],
                                vertices[i],
                                vertices[(i + 1) % len(vertices)],
                                side_long, side_short))
#                                w, h))

  """
  weights[0] = 1 / (weights[3] + weights[0] + weights[1])
  weights[1] = 1 / (weights[0] + weights[1] + weights[2])
  weights[2] = 1 / (weights[1] + weights[2] + weights[3])
  weights[3] = 1 / (weights[2] + weights[3] + weights[0])
  #print "Got residuals", weights
  """

  """
  weights = [0, 0, 0, 0]
  normal = u.cross(other=v).normalize()
  weights[0] = (1 - abs(normal.dot( \
          other=(vertices[0] - c).normalize()))) \
          *    (1 - abs((vertices[3] - vertices[0]).normalize().dot( \
          other=(vertices[1] - vertices[0]).normalize())))
  weights[1] = (1 - abs(normal.dot( \
          other=(vertices[1] - c).normalize()))) \
          *    (1 - abs((vertices[0] - vertices[1]).normalize().dot( \
          other=(vertices[2] - vertices[1]).normalize())))
  weights[2] = (1 - abs(normal.dot(\
          other=(vertices[2] - c).normalize()))) \
          *    (1 - abs((vertices[1] - vertices[2]).normalize().dot( \
          other=(vertices[3] - vertices[2]).normalize())))
  weights[3] = (1 - abs(normal.dot( \
          other=(vertices[3] - c).normalize()))) \
          *    (1 - abs((vertices[2] - vertices[3]).normalize().dot( \
          other=(vertices[0] - vertices[3]).normalize())))
  #print "Got angle-stuff", weights
  """

  #"""
  weights[0] = 1
  weights[1] = 1
  weights[2] = 1
  weights[3] = 1
  #"""

  wtot = weights[0] + weights[1] + weights[2] + weights[3]
  weights[0] /= wtot
  weights[1] /= wtot
  weights[2] /= wtot
  weights[3] /= wtot

  #print "Got weight vector", weights


  """
  t0 = vertices[3] - vertices[0]
  t1 = vertices[1] - vertices[0]
  a  = 180 * t0.angle(t1) / math.pi
  print "GOT %.4f %.4f" % (a, math.fabs(a - 90))

  t0 = vertices[0] - vertices[1]
  t1 = vertices[2] - vertices[1]
  a  = 180 * t0.angle(t1) / math.pi
  print "GOT %.4f %.4f" % (a, math.fabs(a - 90))

  t0 = vertices[1] - vertices[2]
  t1 = vertices[3] - vertices[2]
  a  = 180 * t0.angle(t1) / math.pi
  print "GOT %.4f %.4f" % (a, math.fabs(a - 90))

  t0 = vertices[2] - vertices[3]
  t1 = vertices[0] - vertices[3]
  a  = 180 * t0.angle(t1) / math.pi
  print "GOT %.4f %.4f" % (a, math.fabs(a - 90))
  """

  # THIS IS THE NEW CODE, weighted optimisation, using the unweighted
  # values as starting guesses.  The non-least squares target seems
  # better (more resistant to outliers).  But again tricky and almost
  # useless.
  if False:
    o = optimise_rectangle(vertices, weights, c, u, v, l)
    c = matrix.col((o.x[0], o.x[1], o.x[2]))
    u = matrix.col((o.x[3], o.x[4], o.x[5]))
    v = matrix.col((o.x[6], o.x[7], o.x[8]))
  elif False:
    o = optimise_rectangle_2(vertices, weights, c, u, v, l)
    c = matrix.col((o.x[0], o.x[1], o.x[2]))
    u = matrix.col((o.x[3], o.x[4], o.x[5]))
    v = matrix.col((o.x[6], o.x[7], o.x[8]))
  else:
    pass
    """
    o = optimise_rectangle_3(vertices, matrix.col(weights), u, v, side_long, side_short)
    weights = matrix.col((o.x[0]**2, o.x[1]**2, o.x[2]**2, o.x[3]**2)) / (
      o.x[0]**2 + o.x[1]**2 + o.x[2]**2 + o.x[3]**2)
    u = matrix.col((o.x[4], o.x[5], o.x[6]))
    v = matrix.col((o.x[7], o.x[8], o.x[9]))
    c = (weights[0] * (vertices[0] + v) +
         weights[1] * (vertices[1] - u) +
         weights[2] * (vertices[2] - v) +
         weights[3] * (vertices[3] + u))
    """

  w = (u - v).length()
  h = (u + v).length()
  if w < h:
    t = h
    h = w
    w = t

  vertices_lsq = [c - v, c + u, c + v, c - u]
  rss = weights[0] * (vertices[0] - vertices_lsq[0]).length() \
      + weights[1] * (vertices[1] - vertices_lsq[1]).length() \
      + weights[2] * (vertices[2] - vertices_lsq[2]).length() \
      + weights[3] * (vertices[3] - vertices_lsq[3]).length()
  rss = (vertices[0] - vertices_lsq[0]).norm_sq() \
      + (vertices[1] - vertices_lsq[1]).norm_sq() \
      + (vertices[2] - vertices_lsq[2]).norm_sq() \
      + (vertices[3] - vertices_lsq[3]).norm_sq()


  """
  # This is for printing Gnuplot files for visualisation.

  #print "set arrow from %f,%f,%f to %f,%f,%f nohead" % (tuple(vertices[0]) + tuple(vertices[1]))
  #print "set arrow from %f,%f,%f to %f,%f,%f nohead" % (tuple(vertices[1]) + tuple(vertices[2]))
  #print "set arrow from %f,%f,%f to %f,%f,%f nohead" % (tuple(vertices[2]) + tuple(vertices[3]))
  #print "set arrow from %f,%f,%f to %f,%f,%f nohead" % (tuple(vertices[3]) + tuple(vertices[0]))

  print "%f %f %f" % tuple(vertices[0])
  print "%f %f %f" % tuple(vertices[1])
  print "%f %f %f" % tuple(vertices[2])
  print "%f %f %f" % tuple(vertices[3])
  print "%f %f %f" % tuple(vertices[0])
  print ""
  print "%f %f %f" % tuple(vertices_lsq[0])
  print "%f %f %f" % tuple(vertices_lsq[1])
  print "%f %f %f" % tuple(vertices_lsq[2])
  print "%f %f %f" % tuple(vertices_lsq[3])
  print "%f %f %f" % tuple(vertices_lsq[0])
  """

  """
  print "R: %7.3f %7.3f %7.3f %7.3f -- %d %d -- %f" % (
    (vertices[0] - c + v).length(),
    (vertices[1] - c - u).length(),
    (vertices[2] - c - v).length(),
    (vertices[3] - c + u).length(),
    (u + v).length(), # max((u - v).length(), (u + v).length()),
    (u - v).length(), # min((u - v).length(), (u + v).length()),
    rss)
  print "W: %7.3f %7.3f %7.3f %7.3f -- %f" % (
    weights[0], weights[1], weights[2], weights[3], math.fabs(u.length() - v.length()))
  print "S: %5d %5d %5d %5d" % ((vertices[1] - vertices[0]).length(),
                                (vertices[2] - vertices[1]).length(),
                                (vertices[3] - vertices[2]).length(),
                                (vertices[0] - vertices[3]).length())

  """

  #sys.exit(0)


  return (vertices_lsq, rss)


def _fit_angles(vertices):
  """
  Just a bunch of heuristics.
  """

  # Determine the number of good angles.  An angle is good if it's no
  # more than 0.1 degrees deviation from 90.  Note that it is
  # impossible to have three good angles (because the fourth would
  # have to be good to close the polygon).

  ANGLE_THRESHOLD = 0.01

  good_angles = []
  for i in range(4):
    u = vertices[(i - 1) % 4] - vertices[i]
    v = vertices[(i + 1) % 4] - vertices[i]

    good_angles.append(math.fabs(u.angle(v, deg=True) - 90) < ANGLE_THRESHOLD)
    #print "Vertex", i, "angle", u.angle(v, deg=True), math.fabs(u.angle(v, deg=True) - 90)

  #print "Good angles", good_angles

  if good_angles.count(True) == 0:
    # Must have two parallel edges for this to work (XXX illustrate
    # the non-parallel case).  XXX Need an artificial test case for
    # this!

    a1 = (vertices[1] - vertices[0]).angle(vertices[3] - vertices[2], deg=True)
    a2 = (vertices[3] - vertices[0]).angle(vertices[2] - vertices[1], deg=True)
    if a1 < ANGLE_THRESHOLD:
      print("IN FIX 0a")
      return (vertices, 0)

    elif a2 < ANGLE_THRESHOLD:
      print("IN FIX 0b")
      return (vertices, 0)

  elif good_angles.count(True) == 1:
    # Move the vertex opposite the good angle into its expected
    # position.

    #print "IN FIX 1"

    i = good_angles.index(True)
    u = vertices[(i + 1) % 4] + vertices[(i - 1) % 4] - vertices[i]

    vertices_fixed = copy.deepcopy(vertices)
    vertices_fixed[(i + 2) % 4] = u

    return (vertices_fixed, (vertices[(i + 2) % 4] - u).norm_sq())


  elif good_angles.count(True) == 2:
    # Move the two vertices with bad angles so they line up.  Must be
    # adjacent for this to work (XXX illustrate the non-adjacent case
    # mean)?

    i = good_angles.index(True)
    j = good_angles.index(True, i + 1)
    if i == 0 and j == 3:
      i = 3
      j = 0
    if (i + 1) % 4 != j:
      return (vertices, 0)

    #print "IN FIX 2" # XXX CHECK THIS AGAIN!!

    u = vertices[(i - 1) % 4] - vertices[i]
    v = vertices[(j + 1) % 4] - vertices[j]
    l = 0.5 * (u.length() + v.length())

    u = vertices[i] + l / u.length() * u
    v = vertices[j] + l / v.length() * v

    vertices_fixed = copy.deepcopy(vertices)
    vertices_fixed[(i - 1) % 4] = u
    vertices_fixed[(j + 1) % 4] = v

    return (vertices_fixed,
            (vertices[(i - 1) % 4] - u).norm_sq() +
            (vertices[(j + 1) % 4] - v).norm_sq())


  # Nuttin' to do.
  return (vertices, 0)

def _fit_plane(vertices):
  """
  Fit least-squares plane (see Schomaker et al. (1959) Acta Cryst 12,
  600-604 and Blow (1960) Acta Cryst 13, 168 XXX check WOS).

  @return The residual sum of squares
  """

  # The center is always the mean of the four corners.
  c = 0.25 * (vertices[0] + vertices[1] + vertices[2] + vertices[3])

  # Compute covariance matrix.
  C = matrix.rec(
    elems=tuple(vertices[0] - c) \
      +   tuple(vertices[1] - c) \
      +   tuple(vertices[2] - c) \
      +   tuple(vertices[3] - c),
    n=(4, 3))
  cov = C.transpose_multiply(C)

  # Obtain eigenvalues of covariance matrix.  XXX Ouch, this is ugly.
  from tntbx import svd
  cov_flex = flex.double(flex.grid(3, 3))
  cov_flex += flex.double((cov(0, 0), cov(0, 1), cov(0, 2),
                           cov(1, 0), cov(1, 1), cov(1, 2),
                           cov(2, 0), cov(2, 1), cov(2, 2)))
  svd = svd(cov_flex)

  # Get least-squares residual (smallest eigenvalue) and plane normal
  # (eigenvector corresponding to smallest eigenvalue) from SVD.
  rss = svd.s()[2, 2]
  normal = matrix.col((svd.v()[0, 2], svd.v()[1, 2], svd.v()[2, 2]))

  """
  rss2 = normal.dot(vertices[0] - c)**2 \
      +  normal.dot(vertices[1] - c)**2 \
      +  normal.dot(vertices[2] - c)**2 \
      +  normal.dot(vertices[3] - c)**2
  assert math.fabs(rss - rss2) < 1e-5
  """

  # Project each vertex into the plane.  Normal must have unit length.
  vertices_lsq = [(vertices[0] - c) - normal.dot(vertices[0] - c) * normal + c,
                  (vertices[1] - c) - normal.dot(vertices[1] - c) * normal + c,
                  (vertices[2] - c) - normal.dot(vertices[2] - c) * normal + c,
                  (vertices[3] - c) - normal.dot(vertices[3] - c) * normal + c]

  return (vertices_lsq, rss)



def parse_metrology(path, detector = 'CxiDs1', plot = True, do_diffs = True, old_style_diff_path = None):
  """Measurement file has all units in micrometers.  The coordinate
  system is right-handed XXX verify this.  The origin is at the lower,
  left corner of sensor 1.  XXX Are quadrants still 0, 1, 2, 3 in
  counter-clockwise order beginning at top left corner?"""

  assert detector in ['CxiDs1', 'XppDs1']

  import re

  # Regular expressions for the four kinds of lines that may appear in
  # the metrology definition (XXX measurement file?).  As an
  # extension, everything following a hash mark is ignored.
  # Continuation lines (i.e. backslash preceding a newline) are not
  # permitted.
  pattern_blank = \
      re.compile('^\s*(#.*)?$')
  pattern_quadrant = \
      re.compile('^\s*[Qq][Uu][Aa][Dd]\s+\d+\s*(#.*)?$')
  pattern_sensor = \
      re.compile('^\s*[Ss][Ee][Nn][Ss][Oo][Rr]\s+[Xx]\s+[Yy]\s+[Zz]\s*(#.*)?$')
  pattern_sxyz = \
      re.compile('^\s*\d+\s+[+-]?\d+\s+[+-]?\d+\s+[+-]?\d+\s*(#.*)?$')

  # Parse the file into dict of dict of four-tuples, indexed by
  # quadrant indices and sensor indices, respectively.
  quadrants = {}
  stream = open(path, 'r')
  for line in stream.readlines():
    if pattern_blank.match(line) or pattern_sensor.match(line):
      pass

    elif pattern_quadrant.match(line):
      q = int(line.split()[1])

    elif pattern_sxyz.match(line):
      if 'q' not in locals():
        raise RuntimeError("Unknown quadrant for sensor")

      (s, x, y, z) = [int(i) for i in line.split()[0:4]]
      p = matrix.col((x, y, z))
      if q in quadrants.keys():
        if s in quadrants[q].keys():
          raise RuntimeError("Multiple definitions for sensor")
        quadrants[q][s] = p
      else:
        quadrants[q] = {s: p}

    else:
      raise RuntimeError("Syntax error")
  stream.close()

  if plot:
    def center(coords):
      """ Returns the average of a list of vectors
      @param coords List of vectors to return the center of
      """
      for c in coords:
        if 'avg' not in locals():
          avg = c
        else:
          avg += c
      return avg / len(coords)

    print("Showing un-adjusted, parsed metrology")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    cents = []

    for q_id, quadrant in quadrants.iteritems():
      s = []
      for s_id, sensor_pt in quadrant.iteritems():
        s.append(sensor_pt)
        if len(s) == 4:
          ax.add_patch(Polygon([p[0:2] for p in s], closed=True, color='green', fill=False, hatch='/'))
          cents.append(center(s)[0:2])
          s = []

    for i, c in enumerate(cents):
      ax.annotate("%d:%d"%(i//8,i%8),c)
    ax.set_xlim((0, 200000))
    ax.set_ylim((0, 200000))
    plt.show()

  # More error checking and resolve out-of-order definitions.
  l_diff = []
  s_diff = []
  l_len = []
  s_len = []
  for (q, s) in quadrants.iteritems():
    # Sort the sensor vertices and ensure every sensor has four
    # vertices.
    s_sorted = copy.deepcopy(s.keys())
    s_sorted.sort()
    if len(s_sorted) % 4 != 0:
      raise RuntimeError(
        "Sensor in quadrant does not have integer number of sensors")

    # Ensure the vertices of the sensor are contiguously numbered.
    vertices = []
    for i in range(0, len(s_sorted), 4):
      if    s_sorted[i] + 1 != s_sorted[i + 1] or \
            s_sorted[i] + 2 != s_sorted[i + 2] or \
            s_sorted[i] + 3 != s_sorted[i + 3]:
        raise RuntimeError(
          "Sensor vertices in quadrant not contiguously numbered")
      vertices.append((s[s_sorted[i + 0]],
                       s[s_sorted[i + 1]],
                       s[s_sorted[i + 2]],
                       s[s_sorted[i + 3]]))

    # Organise the quadrant's sensors in XTC stream order.  XXX Open
    # question: should this be done or not?  Should probably require
    # the input to be correctly ordered.  In either case, must ensure
    # first vertex starts a long side.
    if detector == 'CxiDs1': # XXX This applies to the CXI rear detector as well!
      quadrants[q] = _order_sensors(vertices)

    elif detector == 'XppDs1':
      quadrants[q] = {}
      for i in range(len(vertices)):
        v = vertices[i]
        if (v[1] - v[0]).length() < (v[2] - v[1]).length():
          quadrants[q][i] = (v[1], v[2], v[3], v[0])
        else:
          quadrants[q][i] = (v[0], v[1], v[2], v[3])

    # Compute absolute difference between long sides and short sides.
    for i in range(len(vertices)):
      v = vertices[i]

      l = ((v[1] - v[0]).length(), (v[3] - v[2]).length())
      s = ((v[2] - v[1]).length(), (v[0] - v[3]).length())
      if l[0] + l[1] < s[0] + s[1]:
        (l, s) = (s, l)

      l_len.extend(l)
      s_len.extend(s)

      l_diff.append(abs(l[0] - l[1]))
      s_diff.append(abs(s[0] - s[1]))

  l_diff_stat = flex.mean_and_variance(flex.double(l_diff))
  s_diff_stat = flex.mean_and_variance(flex.double(s_diff))
  print("Difference between opposing long  sides: " \
    "%.3f+/-%.3f [%.3f, %.3f] (N = %d)" % (
      l_diff_stat.mean(),
      l_diff_stat.unweighted_sample_standard_deviation(),
      min(l_diff),
      max(l_diff),
      len(l_diff)))
  print("Difference between opposing short sides: " \
    "%.3f+/-%.3f [%.3f, %.3f] (N = %d)" % (
      s_diff_stat.mean(),
      s_diff_stat.unweighted_sample_standard_deviation(),
      min(s_diff),
      max(s_diff),
      len(s_diff)))


  # Get side length distribution parameters with outlier rejection,
  # set SIDE_MAX_SIGMA = float('inf') (in um) to skip iteration.  For
  # stability, stop iterating once no more changes.  XXX What
  # percentage does five sigma correspond to again?
  SIDE_MAX_SIGMA = float('inf')
  (l_mu, l_sigma, l_N, s_mu, s_sigma, s_N) = _find_long_short(quadrants)
  while l_sigma > SIDE_MAX_SIGMA or s_sigma > SIDE_MAX_SIGMA:
    (l_N_old, s_N_old) = (l_N, s_N)
    l_int = (l_mu - 5 * l_sigma, l_mu + 5 * l_sigma)
    s_int = (s_mu - 5 * s_sigma, s_mu + 5 * s_sigma)
    (l_mu, l_sigma, l_N, s_mu, s_sigma, s_N) = _find_long_short(
      quadrants, l_int, s_int)
    if l_N >= l_N_old and s_N >= s_N_old:
      break
  print("Long  side: %.3f+/-%.3f [%.3f, %.3f] (N = %d)" % \
    (l_mu, l_sigma, min(l_len), max(l_len), l_N))
  print("Short side: %.3f+/-%.3f [%.3f, %.3f] (N = %d)" % \
    (s_mu, s_sigma, min(s_len), max(s_len), s_N))

  # XXX Would really expect this to be (2 * 194 + 3) / 185 = 2.11, but
  # it ain't.  Hart et al. (2012) says this should be 20.9 mm by 43.5 mm
  (r_mu, r_sigma, r_N) = _find_aspect_ratio(quadrants)
  print("Aspect ratio: %.3f+/-%.3f (N = %d)" % (r_mu, r_sigma, r_N))


  # Corrections appear to be necessary for 2011-06-20 Ds1 metrology,
  # but no the the 2012-11-08 Dsd ditto.
  corrections = 0

  # For each sensor in each quadrant, determine least-squares fit
  # rectangle with sides (l_mu, s_mu) to the four vertices.
  quadrants_lsq = {}
  for (q, sensors) in quadrants.iteritems():
    print("Q: %1d" % q)
    quadrants_lsq[q] = {}
    for (s, vertices) in sensors.iteritems():

      #if q != 0 or s != 3:
      #  continue

      # XXX This should be optional, really!  It's not really
      # useful (except for quadrant 2, sensor 1).

      vertices_lsq = vertices
      if corrections > 0:
        (vertices_lsq, rss) = _fit_plane(vertices)
        #(vertices_lsq, rss) = _fit_plane(vertices_lsq)
        #print "corrected out-of-plane r.m.s.d.", math.sqrt(0.25 * rss)

      if corrections > 1:
        (vertices_lsq, rss) = _fit_angles(vertices_lsq)
        #print "Corrected bad angles r.m.s.d.", math.sqrt(0.25 * rss)

#      (vertices_lsq, rss) = _fit_rectangle(vertices, l_mu, s_mu)
      (vertices_lsq, rss) = _fit_rectangle(vertices_lsq, l_mu, s_mu)
      quadrants_lsq[q][s] = vertices_lsq
      print("Q: %1d S: %1d r.m.s.d. = %8.3f" % (q, s, math.sqrt(0.25 * rss)))

  # XXX What is aspect ratio after fitting rectangles?
#  (r_mu, r_sigma, r_N) = _find_aspect_ratio(quadrants_lsq)
#  print "Aspect ratio: %.3f+/-%.3f (N = %d)" % (r_mu, r_sigma, r_N)

  # This is a measure of success.
  (l_mu, l_sigma, l_N, s_mu, s_sigma, s_N) = _find_long_short(quadrants_lsq)
  print("Long  side: %.3f+/-%.3f (N = %d)\n" \
      "Short side: %.3f+/-%.3f (N = %d)" % \
      (l_mu, l_sigma, l_N, s_mu, s_sigma, s_N))


  if old_style_diff_path is None:
    return quadrants_lsq

  ##################################
  # METROLOGY REFINEMENT ENDS HERE #
  ##################################


  # The rest is mildly cctbx.xfel specific at the moment.  First,
  # apply transformations: bring to order (slow, fast) <=> (column,
  # row).  This takes care of quadrant rotations, and projects into
  # 2D-space.
  quadrants_trans = {}
  if detector == 'CxiDs1':
    # XXX This applies to the CXI rear detector as well!
    #
    # There is no global origin for the CSPAD:s at CXI, because all
    # four quadrants are measured independently.
    for (q, sensors) in quadrants_lsq.iteritems():
      quadrants_trans[q] = {}

      if q == 0:
        # Q0:
        #   x -> -slow
        #   y -> -fast
        for (s, vertices) in sensors.iteritems():
          quadrants_trans[q][s] = [matrix.col((-v[0], -v[1]))
                                   for v in quadrants_lsq[q][s]]
      elif q == 1:
        # Q1:
        #   x -> +fast
        #   y -> -slow
        for (s, vertices) in sensors.iteritems():
          quadrants_trans[q][s] = [matrix.col((-v[1], +v[0]))
                                   for v in quadrants_lsq[q][s]]
      elif q == 2:
        # Q2:
        #   x -> +slow
        #   y -> +fast
        for (s, vertices) in sensors.iteritems():
          quadrants_trans[q][s] = [matrix.col((+v[0], +v[1]))
                                   for v in quadrants_lsq[q][s]]
      elif q == 3:
        # Q3:
        #   x -> -fast
        #   y -> +slow
        for (s, vertices) in sensors.iteritems():
          quadrants_trans[q][s] = [matrix.col((+v[1], -v[0]))
                                   for v in quadrants_lsq[q][s]]
      else:
        # NOTREACHED
        raise RuntimeError(
          "Detector does not have exactly four quadrants")
  elif detector == 'XppDs1':

    # The XPP CSPAD has fixed quadrants.  For 2013-01-24 measurement,
    # they are all defined in a common coordinate system.
    #
    #   x -> +fast
    #   y -> -slow
    #
    # Fix the global origin to the center of mass of sensor 1 in the
    # four quadrants.
    o = matrix.col((0, 0, 0))
    N = 0
    for (q, sensors) in quadrants_lsq.iteritems():
      for v in sensors[1]:
        o += v
        N += 1
    o /= N

    for (q, sensors) in quadrants_lsq.iteritems():
      quadrants_trans[q] = {}
      for (s, vertices) in sensors.iteritems():
        quadrants_trans[q][s] = [matrix.col((-v[1] - (-o[1]), +v[0] - (+o[0])))
                                 for v in quadrants_lsq[q][s]]
  else:
    raise RuntimeError(
      "Unknown convention '%s'" % detector)

  # Read old-style metrology.
  from xfel.cxi.cspad_ana.parse_calib import calib2sections
  calib_dir = old_style_diff_path

  sections = calib2sections(calib_dir)

  aa_new = []

  # Now working with (slow, fast) in 2D.  Enforce right angles and
  # integer pixel positions.  Temporary thingy: output the vertices in
  # Q0 convention to stdout.  Scale factor (micrometers per pixel).
  if plot:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')

  pixel_size = 109.92 # Correct pixel size at 20 degrees (Hart et al., 2012).
  corrections = {}

  # XXX Whoa!!  Either (1, 1) or (1, -1), possibly (-1, 1).
  sign_hack_s = (1, 1)
  sign_hack_l = sign_hack_s

  for (q, sensors) in quadrants_trans.iteritems():
    corrections[q] = {}

    # Figure out quadrant offset from old-style data.  Get bottom ASIC
    # of sensor 1.  Lower right corner of all the unrotated quadrants
    # is always origin.  XXX This is really a bunch of
    # magic-number-hacks for February 2013 Dsd to make things look
    # nice--we don't know exactly where the ASIC corner is with
    # respect to the nearest measurement point.
    offset = (0, 0)
    c = sections[q][1].corners_asic()[0]

    if detector == 'CxiDs1':
      if q == 0:
        offset = (c[2] + 2, c[3] + 3)
      elif q == 1:
        offset = (c[2] + 3, c[1] - 3)
      elif q == 2:
        offset = (c[0] - 3, c[1] - 2)
      elif q == 3:
        offset = (c[0] - 3, c[3] + 3)

    for (s, vertices) in sensors.iteritems():
      corrections[q][s] = {}

      # Compute differences against old-style metrology.
      corners_old = sections[q][s].corners_asic()

      # l is the unit vector oriented along the sensor's long side
      # (from the first ASIC to the second).  ns will be the unit
      # vector oriented along the sensor's short side.  XXX Should
      # this reordering have been resolved above?  The sensor rotation
      # is *very* poorly defined at the moment!  Sensors 6 and 7 do
      # not follow the pattern!
#      if s == 6 or s == 7:
#        l = (vertices[0] - vertices[1]).normalize()
#      else:
#        l = (vertices[1] - vertices[0]).normalize()

      l = (vertices[1] - vertices[0]).normalize()
      if   q == 0 and s in [0, 1]       or \
           q == 2 and s in [4, 5]       or \
           q == 3 and s in [2, 3, 6, 7]:
        # l should point in -slow direction
        if l.elems[0] > 0:
          l *= -1
      elif q == 0 and s in [2, 3, 6, 7] or \
           q == 1 and s in [0, 1]       or \
           q == 3 and s in [4, 5]:
        # l should point in +fast direction
        if l.elems[1] < 0:
          l *= -1
      elif q == 0 and s in [4, 5]       or \
           q == 1 and s in [2, 3, 6, 7] or \
           q == 2 and s in [0, 1]:
        # l should point in +slow direction
        if l.elems[0] < 0:
          l *= -1
      elif q == 1 and s in [4, 5]       or \
           q == 2 and s in [2, 3, 6, 7] or \
           q == 3 and s in [0, 1]:
        # l should point in -fast direction
        if l.elems[1] > 0:
          l *= -1
      else:
        # NOTREACHED
        raise RuntimeError(
          "Unknown direction for sensor %d in quadrant %d" % (s, q))

      # Force l to right angle
      if l.elems[0] < -0.5:
        l = matrix.col((-1, 0))
      elif l.elems[0] > +0.5:
        l = matrix.col((+1, 0))
      elif l.elems[1] < -0.5:
        l = matrix.col((0, -1))
      elif l.elems[1] > +0.5:
        l = matrix.col((0, +1))
      else:
        print("This should not happen!")

      # According to Henrik Lemke, the XPP detector is actually
      # rotated by 180 degrees with respect to the optical metrology
      # measurements.
      if detector == 'XppDs1':
        l *= -1
        for i in range(len(vertices)):
          vertices[i] *= -1

#      print "Quadrant %d sensor %d vector l [% 3.2f, % 3.2f]" % \
#        (q, s, l.elems[0], l.elems[1])

      # m is the 2D center of mass, in units of pixels.  XXX With the
      # most recent change, c is no longer used?  Done to ensure
      # ASIC:s on sensor align, and three-pixel gap is correct (even
      # though rotations are completely discarded).
      m = 0.25 * (vertices[0] +
                  vertices[1] +
                  vertices[2] +
                  vertices[3]) / pixel_size

      if abs(l.elems[0]) > abs(l.elems[1]):
        # Sensor is standing up.

        # First ASIC, now in units of pixels.  XXX Note, no longer
        # considering center of pixels.  Pixel coordinate is at top
        # left corner of pixel?  XXX Note: this aint't always the top
        # or the bottom ASIC.
        c = m - (194 + 3) / 2 * l
        corners_new = [
          int(math.floor(c.elems[0] - (194 - 1) / 2)) + offset[0],
          int(math.floor(c.elems[1] - (185 - 1) / 2)) + offset[1],
          int(math.ceil( c.elems[0] + (194 - 1) / 2)) + offset[0],
          int(math.ceil( c.elems[1] + (185 - 1) / 2)) + offset[1]]
        assert  corners_new[2] - corners_new[0] == 194 \
            and corners_new[3] - corners_new[1] == 185

        # XXX Sign of correction!
        corrections_tl = [sign_hack_s[0] * (corners_new[0] - corners_old[0][0]),
                          sign_hack_s[1] * (corners_new[1] - corners_old[0][1])]
        corrections_br = [sign_hack_s[0] * (corners_new[2] - corners_old[0][2]),
                          sign_hack_s[1] * (corners_new[3] - corners_old[0][3])]
        assert corrections_tl == corrections_br # XXX for test
        corrections[q][s][0] = corrections_tl

        if q == 0 and s == 0:
          print("HATTNE diag PRUTT #0", corrections_tl)

        if plot:
          v_new = ((corners_new[0], corners_new[1]),
                   (corners_new[2], corners_new[1]),
                   (corners_new[2], corners_new[3]),
                   (corners_new[0], corners_new[3]))
          v_old = ((corners_old[0][0], corners_old[0][1]),
                   (corners_old[0][2], corners_old[0][1]),
                   (corners_old[0][2], corners_old[0][3]),
                   (corners_old[0][0], corners_old[0][3]))
          ax.add_patch(Polygon(
              v_new, closed=True, color='green', fill=False, hatch='/'))
          ax.add_patch(Polygon(
              v_old, closed=True, color='red', fill=False))

          aa_new += [corners_new[0], corners_new[1],
                     corners_new[2], corners_new[3]]

        # Second ASIC, now in units of pixels.  XXX Note: this aint't
        # always the top or the bottom ASIC.

        c = m + (194 + 3) / 2 * l
        corners_new = [
          int(math.floor(c.elems[0] - (194 - 1) / 2)) + offset[0],
          int(math.floor(c.elems[1] - (185 - 1) / 2)) + offset[1],
          int(math.ceil( c.elems[0] + (194 - 1) / 2)) + offset[0],
          int(math.ceil( c.elems[1] + (185 - 1) / 2)) + offset[1]]
        assert  corners_new[2] - corners_new[0] == 194 \
            and corners_new[3] - corners_new[1] == 185

        # XXX Sign of correction!
        corrections_tl = [sign_hack_s[0] * (corners_new[0] - corners_old[1][0]),
                          sign_hack_s[1] * (corners_new[1] - corners_old[1][1])]
        corrections_br = [sign_hack_s[0] * (corners_new[2] - corners_old[1][2]),
                          sign_hack_s[1] * (corners_new[3] - corners_old[1][3])]
        assert corrections_tl == corrections_br # XXX for test
        corrections[q][s][1] = corrections_tl

        if q == 0 and s == 0:
          print("HATTNE diag PRUTT #1", corrections_tl)


        if plot:
          v_new = ((corners_new[0], corners_new[1]),
                   (corners_new[2], corners_new[1]),
                   (corners_new[2], corners_new[3]),
                   (corners_new[0], corners_new[3]))
          v_old = ((corners_old[1][0], corners_old[1][1]),
                   (corners_old[1][2], corners_old[1][1]),
                   (corners_old[1][2], corners_old[1][3]),
                   (corners_old[1][0], corners_old[1][3]))
          ax.add_patch(Polygon(
              v_new, closed=True, color='green', fill=False, hatch='/'))
          ax.add_patch(Polygon(
              v_old, closed=True, color='red', fill=False))

          aa_new += [corners_new[0], corners_new[1],
                     corners_new[2], corners_new[3]]

      else:
        # Sensor is laying down.

        # First ASIC (not necessarily the left), now in units of
        # pixels .
        c = m - (194 + 3) / 2 * l
        corners_new = [
          int(math.floor(c.elems[0] - (185 - 1) / 2)) + offset[0],
          int(math.floor(c.elems[1] - (194 - 1) / 2)) + offset[1],
          int(math.ceil( c.elems[0] + (185 - 1) / 2)) + offset[0],
          int(math.ceil( c.elems[1] + (194 - 1) / 2)) + offset[1]]
        assert  corners_new[2] - corners_new[0] == 185 \
            and corners_new[3] - corners_new[1] == 194

        # XXX Sign of correction!
        corrections_tl = [sign_hack_l[0] * (corners_new[0] - corners_old[0][0]),
                          sign_hack_l[1] * (corners_new[1] - corners_old[0][1])]
        corrections_br = [sign_hack_l[0] * (corners_new[2] - corners_old[0][2]),
                          sign_hack_l[1] * (corners_new[3] - corners_old[0][3])]
        assert corrections_tl == corrections_br # XXX for test
        corrections[q][s][0] = corrections_tl

#        if q == 0 and s == 6:
#          print "HATTNE diag #0", corrections_tl

        if plot:
          v_new = ((corners_new[0], corners_new[1]),
                   (corners_new[2], corners_new[1]),
                   (corners_new[2], corners_new[3]),
                   (corners_new[0], corners_new[3]))
          v_old = ((corners_old[0][0], corners_old[0][1]),
                   (corners_old[0][2], corners_old[0][1]),
                   (corners_old[0][2], corners_old[0][3]),
                   (corners_old[0][0], corners_old[0][3]))
          ax.add_patch(Polygon(
              v_new, closed=True, color='green', fill=False, hatch='/'))
          ax.add_patch(Polygon(
              v_old, closed=True, color='red', fill=False))

          aa_new += [corners_new[0], corners_new[1],
                     corners_new[2], corners_new[3]]

        # Second ASIC (not necessarily the right), now in units of
        # pixels.
        c = m + (194 + 3) / 2 * l
        corners_new = [
          int(math.floor(c.elems[0] - (185 - 1) / 2)) + offset[0],
          int(math.floor(c.elems[1] - (194 - 1) / 2)) + offset[1],
          int(math.ceil( c.elems[0] + (185 - 1) / 2)) + offset[0],
          int(math.ceil( c.elems[1] + (194 - 1) / 2)) + offset[1]]
        assert  corners_new[2] - corners_new[0] == 185 \
            and corners_new[3] - corners_new[1] == 194

        # XXX Sign of correction!
        corrections_tl = [sign_hack_l[0] * (corners_new[0] - corners_old[1][0]),
                          sign_hack_l[1] * (corners_new[1] - corners_old[1][1])]
        corrections_br = [sign_hack_l[0] * (corners_new[2] - corners_old[1][2]),
                          sign_hack_l[1] * (corners_new[3] - corners_old[1][3])]
        assert corrections_tl == corrections_br # XXX for test
        corrections[q][s][1] = corrections_tl

#        if q == 0 and s == 6:
#          print "HATTNE diag #1", corrections_tl

        if plot:
          v_new = ((corners_new[0], corners_new[1]),
                   (corners_new[2], corners_new[1]),
                   (corners_new[2], corners_new[3]),
                   (corners_new[0], corners_new[3]))
          v_old = ((corners_old[1][0], corners_old[1][1]),
                   (corners_old[1][2], corners_old[1][1]),
                   (corners_old[1][2], corners_old[1][3]),
                   (corners_old[1][0], corners_old[1][3]))
          ax.add_patch(Polygon(
              v_new, closed=True, color='green', fill=False, hatch='/'))
          ax.add_patch(Polygon(
              v_old, closed=True, color='red', fill=False))

          aa_new += [corners_new[0], corners_new[1],
                     corners_new[2], corners_new[3]]

  if plot:
    #ax.set_xlim((0, 1850))
    #ax.set_ylim((0, 1850))
    ax.set_xlim((-1000, 2000))
    ax.set_ylim((-1000, 2000))
    plt.show()

  # Build list and output differences to old-style metrology.
  corrections_list = flex.int(4 * 8 * 2 * 2)
  for (q, sensors) in corrections.iteritems():
    for (s, correction) in sensors.iteritems():
      corrections_list[q * 8 * 2 * 2 + s * 2 * 2 + 0 * 2 + 0] = correction[0][0]
      corrections_list[q * 8 * 2 * 2 + s * 2 * 2 + 0 * 2 + 1] = correction[0][1]
      corrections_list[q * 8 * 2 * 2 + s * 2 * 2 + 1 * 2 + 0] = correction[1][0]
      corrections_list[q * 8 * 2 * 2 + s * 2 * 2 + 1 * 2 + 1] = correction[1][1]
  print("corrected_auxiliary_translations =", list(corrections_list))

  for i in range(len(aa_new)):
    if detector == 'XppDs1':
      # XXX Not actually tested for XPP.  XXX magic number, again!
      aa_new[i] += 1765 // 2

    assert aa_new[i] >= 0 and aa_new[i] <= 1765

  print("new active areas", len(aa_new), aa_new)


  if plot:
    print("Showing active areas")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    cents = []

    for item in [(aa_new[(i*4)+0],
                 aa_new[(i*4)+1],
                 aa_new[(i*4)+2],
                 aa_new[(i*4)+3]) for i in range(64)]:
      x1,y1,x2,y2 = item
      tile = [(x1,y1),(x2,y1),(x2,y2),(x1,y2)]
      ax.add_patch(Polygon(tile, closed=True, color='green', fill=False, hatch='/'))
      cents.append(center([matrix.col(p) for p in tile])[0:2])

    for i, c in enumerate(cents):
      ax.annotate(i,c)
    ax.set_xlim((0, 2000))
    ax.set_ylim((0, 2000))
    plt.show()


  return quadrants_lsq


# phenix.python flatfile.py 2011-08-10-Metrology.txt
if __name__ == '__main__':
  import libtbx.load_env
  assert len(sys.argv[1:]) == 1
  parse_metrology(sys.argv[1],
                  old_style_diff_path=libtbx.env.find_in_repositories("xfel/metrology/CSPad/run4/CxiDs1.0_Cspad.0"))
