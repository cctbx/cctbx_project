from __future__ import division
import smtbx.refinement.constraints as _
from smtbx.refinement.constraints import InvalidConstraint
import math
from scitbx.math import superpose
from scitbx.array_family import flex

class rigid_pivoted_rotable_group(object):
  """ a set of atoms (rigid body) rides on a pivot atom and rotates around
  the pivot-pivot_neighbour bond, the original geometry is not altered
  """

  def __init__(self, pivot, pivot_neighbour, ind_sequence, sizeable, rotable):
    if len(ind_sequence) == 0:
      raise InvalidConstraint("at least one atom is expected")
    self.pivot = pivot
    self.pivot_neighbour = pivot_neighbour
    self.indices = ind_sequence
    self.sizeable = sizeable
    self.rotable = rotable

  def __eq__(self, other):
    if (self.pivot != other.pivot or\
        self.pivot_neighbour != other.pivot_neghbour or\
        self.indices != other.indices):
      return False
    return True

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    pivot_sp = reparametrisation.add_new_site_parameter(self.pivot)
    pivot_n_sp = reparametrisation.add_new_site_parameter(self.pivot_neighbour)
    azimuth = reparametrisation.add(_.independent_scalar_parameter,
                                    value=0, variable=self.rotable)
    size = reparametrisation.add(_.independent_scalar_parameter,
                                    value=1, variable=self.sizeable)
    scatterers = tuple([scatterers[i] for i in self.indices])
    param = reparametrisation.add(
      _.rigid_pivoted_rotable_group,
      pivot=pivot_sp,
      pivot_neighbour=pivot_n_sp,
      azimuth=azimuth,
      size=size,
      scatterers=scatterers)
    for i, j in enumerate(self.indices):
      reparametrisation.add_new_site_proxy_parameter(param, i, j)
      reparametrisation.asu_scatterer_parameters[j].site = param

class rigid_rotable_expandable_group(object):
  """ a set of atoms rides on a pivot atom, rotates and uniformly
  expands or shrinks
  """

  def __init__(self, center, ind_sequence, sizeable, rotable):
    if len(ind_sequence) == 0:
      raise InvalidConstraint("at least one atom is expected")
    self.pivot = center
    self.indices = ind_sequence
    self.sizeable = sizeable
    self.rotable = rotable

  def __eq__(self, other):
    if (self.pivot != other.pivot or self.indices != other.indices):
      return False
    return True

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    pivot_sp = reparametrisation.add_new_site_parameter(self.pivot)
    size = reparametrisation.add(_.independent_scalar_parameter,
                                    value=1, variable=self.sizeable)
    alpha = reparametrisation.add(_.independent_scalar_parameter,
                                    value=0, variable=self.rotable)
    beta = reparametrisation.add(_.independent_scalar_parameter,
                                    value=0, variable=self.rotable)
    gamma = reparametrisation.add(_.independent_scalar_parameter,
                                    value=0, variable=self.rotable)
    scatterers = tuple([scatterers[i] for i in self.indices])
    param = reparametrisation.add(
      _.rigid_rotable_expandable_group,
      pivot=pivot_sp,
      size=size,
      alpha=alpha,
      beta=beta,
      gamma=gamma,
      scatterers = scatterers)
    for i, j in enumerate(self.indices):
      reparametrisation.add_new_site_proxy_parameter(param, i, j)
      reparametrisation.asu_scatterer_parameters[j].site = param


class idealised_fragment(object):
  """ ported from olex2 xlib/fragment.h
  generates parameterised coordinates for four framents:
  Cp, Ph, Cp* and naphthalene
  """
  def __init__(self):
    self.default_lengths = {
      "Cp"  : (1.42,),
      "Cp*" : (1.42, 1.063),
      "Ph"  : (1.39,),
      "Naphthalene" : (1.39,)
      }
  class point:  #helper class
    def __init__(self, x, y):
      self.x = x
      self.y = y
    def __mul__(self, k):
      return idealised_fragment.point(self.x*k, self.y*k)
    def __add__(self, p):
      return idealised_fragment.point(self.x+p.x, self.y+p.y)
    def length(self):
      return math.sqrt(self.x*self.x+self.y*self.y)
    def __repr__(self):
      return "(" + str(self.x) + "," + str(self.y) + ")"

  def generate_ring(self, edges, r):
    """ generates coordinates of a ring with given radius
    """
    angle = math.pi*2/edges
    ca = math.cos(angle)
    sa = math.sin(angle)
    result = []
    p = idealised_fragment.point(ca, -sa)
    for i in xrange(0, edges):
      result.append(p*r)
      x = p.x
      p.x = ca*x + sa*p.y
      p.y = ca*p.y - sa*x
    return result

  def generate_fragment(self, fragment, lengths=None):
    """ generates given fragment with given/default bond lengths
    returns a list of idealised_fragment.point, having x and y
    attributes
    """
    if lengths == None:
      lengths = self.default_lengths[fragment]
    if fragment == "Cp":
      return self.generate_ring(5, 0.5*lengths[0]/math.cos(54*math.pi/180))
    if fragment == "Ph":
      return self.generate_ring(6, lengths[0])
    if fragment == "Cp*":
      r = 0.5*lengths[0]/math.cos(54*math.pi/180)
      res = self.generate_ring(5, r)
      for i in self.generate_ring(5, r+lengths[1]):
        res.append(i)
      return res
    if fragment == "Naphthalene":
      res = self.generate_ring(6, lengths[0])
      res.append(res[0])
      for i in xrange(3,6):  res.append(res[i]);
      p = res[4]+res[5]
      p = p * (lengths[0]*2*math.cos(math.pi/6)/p.length())
      for i in xrange(6, len(res)):
        res[i] = res[i] + p
      p = res[7]
      res[7] = res[9]
      res[9] = p
      return res

  def fit(self, fragment, reference_sites, control_point_indices=None):
    """ fits given fragment to given sites, if control_points indices
    are not given - all points are fit, otherwise only control points
    are fit and the result is propagated to the rest of the fragment
    coordinates. returns coordinates of the trasformed fragment
    """
    if not control_point_indices:
      control_point_indices = range(0, len(fragment))
    to_fit = [
      (fragment[i].x, fragment[i].y, 0) for i in control_point_indices]
    lsf = superpose.least_squares_fit(
      flex.vec3_double(reference_sites), flex.vec3_double(to_fit))
    to_fit = flex.vec3_double([(i.x, i.y, 0) for i in fragment])
    return lsf.r.elems * to_fit + lsf.t.elems
