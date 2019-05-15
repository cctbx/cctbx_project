from __future__ import absolute_import, division, print_function
import smtbx.refinement.constraints as _
from smtbx.refinement.constraints import InvalidConstraint
import math
from scitbx.math import superpose
from scitbx import matrix
from scitbx.array_family import flex
import itertools
from six.moves import range

class rigid_pivoted_rotatable_group(object):
  """ a set of atoms (rigid body) rides on a pivot atom and rotates around
  the pivot-pivot_neighbour bond, the original geometry is not altered
  """

  def __init__(self, pivot, pivot_neighbour, ind_sequence, sizeable, rotatable):
    if len(ind_sequence) == 0:
      raise InvalidConstraint("at least one atom is expected")
    self.pivot = pivot
    self.pivot_neighbour = pivot_neighbour
    self.indices = ind_sequence
    self.sizeable = sizeable
    self.rotatable = rotatable

  @property
  def constrained_parameters(self):
    return tuple((idx, 'site') for idx in self.indices)

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    pivot_sp = reparametrisation.add_new_site_parameter(self.pivot)
    pivot_n_sp = reparametrisation.add_new_site_parameter(self.pivot_neighbour)
    azimuth = reparametrisation.add(_.independent_scalar_parameter,
                                    value=0, variable=self.rotatable)
    size = reparametrisation.add(_.independent_scalar_parameter,
                                    value=1, variable=self.sizeable)
    scatterers = tuple([scatterers[i] for i in self.indices])
    param = reparametrisation.add(
      _.rigid_pivoted_rotatable_group,
      pivot=pivot_sp,
      pivot_neighbour=pivot_n_sp,
      azimuth=azimuth,
      size=size,
      scatterers=scatterers)
    for i, j in enumerate(self.indices):
      reparametrisation.add_new_site_proxy_parameter(param, i, j)
      reparametrisation.asu_scatterer_parameters[j].site = param

class rigid_rotatable_expandable_group(object):
  """ a set of atoms rides on a pivot atom, rotates and uniformly
  expands or contracts
  """

  def __init__(self, center, ind_sequence, sizeable, rotatable):
    if len(ind_sequence) == 0:
      raise InvalidConstraint("at least one atom is expected")
    self.pivot = center
    self.indices = ind_sequence
    self.sizeable = sizeable
    self.rotatable = rotatable

  @property
  def constrained_parameters(self):
    return tuple((idx, 'site') for idx in self.indices)

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    pivot_sp = reparametrisation.add_new_site_parameter(self.pivot)
    size = reparametrisation.add(_.independent_scalar_parameter,
                                    value=1, variable=self.sizeable)
    alpha = reparametrisation.add(_.independent_scalar_parameter,
                                    value=0, variable=self.rotatable)
    beta = reparametrisation.add(_.independent_scalar_parameter,
                                    value=0, variable=self.rotatable)
    gamma = reparametrisation.add(_.independent_scalar_parameter,
                                    value=0, variable=self.rotatable)
    scatterers = tuple([scatterers[i] for i in self.indices])
    param = reparametrisation.add(
      _.rigid_rotatable_expandable_group,
      pivot=pivot_sp,
      size=size,
      alpha=alpha,
      beta=beta,
      gamma=gamma,
      scatterers = scatterers)
    for i, j in enumerate(self.indices):
      reparametrisation.add_new_site_proxy_parameter(param, i, j)
      reparametrisation.asu_scatterer_parameters[j].site = param


class rigid_riding_expandable_group(object):
  """ a set of atoms rides on a pivot atom, rotates and uniformly
  expands or shrinks
  """

  def __init__(self, center, ind_sequence, sizeable):
    if len(ind_sequence) == 0:
      raise InvalidConstraint("at least one atom is expected")
    self.pivot = center
    self.indices = ind_sequence
    self.sizeable = sizeable

  @property
  def constrained_parameters(self):
    return tuple((idx, 'site') for idx in self.indices)

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    pivot_sp = reparametrisation.add_new_site_parameter(self.pivot)
    size = reparametrisation.add(_.independent_scalar_parameter,
                                    value=1, variable=self.sizeable)
    scatterers = tuple([scatterers[i] for i in self.indices])
    param = reparametrisation.add(
      _.rigid_riding_expandable_group,
      pivot=pivot_sp,
      size=size,
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
    for i in range(0, edges):
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
      for i in range(3,6):  res.append(res[i]);
      p = res[4]+res[5]
      p = p * (lengths[0]*2*math.cos(math.pi/6)/p.length())
      for i in range(6, len(res)):
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


class same_group(object):
  """ non-crystallographic symmetry constraint
  """

  def __init__(self, groups, angles=(0,0,0), fix_xyz=True, fix_u=True):
    """ fix_xyz and fix_u are to be used for debugging purposes only:
    if the coordinates are not fixed, then the refined angles should
    be stored externally and passed to this object
    """
    if len(groups) < 2:
      raise InvalidConstraint("at least two atom sets are expected")
    l = len(groups[0])
    for g in groups[1:]:
      if len(g) != l:
        raise InvalidConstraint("atoms sets differ in size")
    self.groups = groups
    self.fix_xyz = fix_xyz
    self.fix_u = fix_u
    self.angles = angles

  @property
  def constrained_parameters(self):
    result = ()
    for g in itertools.islice(self.groups, 1):
      for i in g:
        if self.fix_xyz:
          result += ((i, 'site'),)
        if self.fix_u:
          result += ((i, 'U'),)
    return result

  def add_to(self, reparametrisation):
    if not self.fix_u and not self.fix_xyz:
      return
    scatterers = reparametrisation.structure.scatterers()
    ref_sites = []
    ref_u_isos = []
    ref_u_stars = []
    ref_adps = []
    src_crds = []
    inv_src_crds = []
    uc = reparametrisation.structure.unit_cell()
    for i in self.groups[0]:
      src_crds.append(uc.orthogonalize(scatterers[i].site))
      if self.fix_xyz:
        ref_sites.append(reparametrisation.add_new_site_parameter(i))
      if self.fix_u:
        if scatterers[i].flags.use_u_iso():
          ref_u_isos.append(
            reparametrisation.add_new_thermal_displacement_parameter(i))
        else:
          ref_u_stars.append(
            reparametrisation.add_new_thermal_displacement_parameter(i))
    for g in self.groups[1:]:
      if len(g) != len(self.groups[0]):
        raise InvalidConstraint("Group size mismatch")
      g_scatterers = []
      g_u_iso_scatterers  =[]
      g_u_star_scatterers = []
      crds = []
      for idx, i in enumerate(g):
        if scatterers[i].flags.use_u_iso() !=\
           scatterers[self.groups[0][idx]].flags.use_u_iso():
          raise InvalidConstraint("Mixing isotropic and anisotropic parameters")
        g_scatterers.append(scatterers[i])
        crds.append(uc.orthogonalize(scatterers[i].site))
        if scatterers[i].flags.use_u_iso():
          g_u_iso_scatterers.append(scatterers[i])
        else:
          g_u_star_scatterers.append(scatterers[i])
      #need to map reference to target
      lsf = superpose.least_squares_fit(
        flex.vec3_double(crds), flex.vec3_double(src_crds))
      #create a list of inverted coordinates if needed
      if len(inv_src_crds) == 0:
        for i in range(0, len(g)):
          inv_src_crds.append(
            2*matrix.col(lsf.other_shift)-matrix.col(src_crds[i]))
      rm = lsf.r
      t = matrix.col(lsf.reference_shift)-matrix.col(lsf.other_shift)
      new_crd = lsf.other_sites_best_fit()
      d = 0
      for i, c in enumerate(new_crd):
        d += matrix.col(matrix.col(c)-matrix.col(crds[i])).length_sq()
      lsf = superpose.least_squares_fit(
        flex.vec3_double(crds), flex.vec3_double(inv_src_crds))
      new_crd = lsf.other_sites_best_fit()
      d_inv = 0
      for i, c in enumerate(new_crd):
        d_inv += matrix.col(matrix.col(c)-matrix.col(crds[i])).length_sq()
      if d_inv < d:
        rm = -lsf.r
      if self.fix_xyz:
        shifts_and_angles =\
          reparametrisation.add(_.independent_small_6_vector_parameter,
                                value=(t[0],t[1],t[2],0,0,0), variable=True)
        if len(ref_u_stars) > 0:
          u_star_param = reparametrisation.add(
            _.same_group_u_star,
            scatterers=g_u_star_scatterers,
            u_stars=ref_u_stars,
            alignment_matrix=rm,
            shifts_and_angles=shifts_and_angles
          )
      elif len(ref_u_stars) > 0:
        angles =\
          reparametrisation.add(_.independent_small_3_vector_parameter,
                                value=self.angles, variable=True)
        u_star_param = reparametrisation.add(
          _.same_group_u_star,
          scatterers=g_u_star_scatterers,
          u_stars=ref_u_stars,
          alignment_matrix=rm,
          angles=angles
        )
      if self.fix_xyz:
        site_param = reparametrisation.add(
          _.same_group_xyz,
          scatterers=g_scatterers,
          sites=ref_sites,
          alignment_matrix=rm,
          shifts_and_angles=shifts_and_angles
        )
      if len(ref_u_isos) > 0:
        u_iso_param = reparametrisation.add(
          _.same_group_u_iso,
          scatterers=g_u_iso_scatterers,
          u_isos=ref_u_isos
        )
      site_proxy_index = 0
      u_star_proxy_index = 0
      u_iso_proxy_index = 0
      for i in g:
        if self.fix_xyz:
          reparametrisation.asu_scatterer_parameters[i].site = site_param
          reparametrisation.add_new_same_group_site_proxy_parameter(
            site_param, site_proxy_index, i)
          site_proxy_index += 1
        if self.fix_u:
          if scatterers[i].flags.use_u_iso():
            reparametrisation.asu_scatterer_parameters[i].u = u_iso_param
            reparametrisation.shared_Us[i] = reparametrisation.add(
              _.same_group_u_iso_proxy,
              parent=u_iso_param,
              index=u_iso_proxy_index
              )
            u_iso_proxy_index += 1
          else:
            reparametrisation.asu_scatterer_parameters[i].u = u_star_param
            reparametrisation.shared_Us[i] = reparametrisation.add(
              _.same_group_u_star_proxy,
              parent=u_star_param,
              index=u_star_proxy_index
              )
            u_star_proxy_index += 1
