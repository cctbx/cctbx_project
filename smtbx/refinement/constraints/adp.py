from __future__ import absolute_import, division, print_function
import smtbx.refinement.constraints as _
from smtbx.refinement.constraints import InvalidConstraint
from math import pi
from six.moves import range

class u_iso_proportional_to_pivot_u_eq(object):
  """ u_iso of some scatterer constrained to be proportional to
      equivalent u_iso associated with adp of another scatterer
  """

  __slots__ = ('u_iso_scatterer_idx', 'u_eq_scatterer_idx', 'multiplier')

  def __init__(self, u_iso_scatterer_idx, u_eq_scatterer_idx, multiplier):
    self.u_iso_scatterer_idx = u_iso_scatterer_idx
    self.u_eq_scatterer_idx = u_eq_scatterer_idx
    self.multiplier = multiplier

  def __eq__(self, other):
    """ For debugging purposes mostly as it is not needed by the
        constraint framework.
    """
    return (self.u_iso_scatterer_idx == other.u_iso_scatterer_idx
            and self.u_eq_scatterer_idx == other.u_eq_scatterer_idx
            and self.multiplier == other.multiplier)

  @property
  def constrained_parameters(self):
    return tuple(((self.u_iso_scatterer_idx, 'U'),))

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    if scatterers[self.u_eq_scatterer_idx].flags.use_u_aniso():
      param = reparametrisation.add(
        _.u_iso_proportional_to_pivot_u_eq,
        pivot_u=reparametrisation.add_new_thermal_displacement_parameter(
          self.u_eq_scatterer_idx),
        scatterer = reparametrisation.structure.scatterers()[
          self.u_iso_scatterer_idx],
        multiplier=self.multiplier)
    else:
      param = reparametrisation.add(
        _.u_iso_proportional_to_pivot_u_iso,
        pivot_u_iso=reparametrisation.add_new_thermal_displacement_parameter(
          self.u_eq_scatterer_idx),
        scatterer = scatterers[self.u_iso_scatterer_idx],
        multiplier=self.multiplier)
    reparametrisation.asu_scatterer_parameters[
      self.u_iso_scatterer_idx].u = param

class shared_u(object):
  """ u_iso or u_star of some scatterer constrained to be equal to
      u_iso or u_star of another scatterer
  """

  def __init__(self, ind_sequence):
    if len(ind_sequence) < 2:
      raise InvalidConstraint("at least two atoms are expected")
    self.indices = ind_sequence

  @property
  def constrained_parameters(self):
    return tuple((idx, "U") for idx in self.indices[1:])

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    src_uses_u = scatterers[self.indices[0]].flags.use_u_aniso()
    for i in range(1, len(self.indices)):
      if scatterers[self.indices[i]].flags.use_u_aniso() != src_uses_u:
        raise InvalidConstraint(
          "mixing isotropic and anisotropic atoms is not allowed for shared ADP")

    u_c = reparametrisation.add_new_thermal_displacement_parameter(
      self.indices[0])
    for i in range(1, len(self.indices)):
      if src_uses_u:
        param = reparametrisation.add(
          _.shared_u_star,
          reference=u_c,
          scatterer = scatterers[self.indices[i]])
      else:
        param = reparametrisation.add(
          _.shared_u_iso,
          reference=u_c,
          scatterer = scatterers[self.indices[i]])
      reparametrisation.shared_Us[self.indices[i]] = u_c
      reparametrisation.asu_scatterer_parameters[self.indices[i]].u = param
    self.value = u_c

class shared_rotated_u(object):
  """ u_eq or u_star of some scatterer constrained to be equal to
      u_iso or u_start of another scatterer
  """

  def __init__(self, ind_ref, ind_atom, direction,
               angle_value, refine_angle=False):
    self.ind_ref = ind_ref
    self.ind_atom = ind_atom
    self.direction = direction
    self.angle_value = angle_value
    self.refine_angle = bool(refine_angle)

  @property
  def constrained_parameters(self):
    return tuple(((self.ind_atom, 'U'),))

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    if not scatterers[self.ind_ref].flags.use_u_aniso() or\
       not scatterers[self.ind_atom].flags.use_u_aniso():
      raise InvalidConstraint(
        "only anisotropic atoms are allowed for shared rotated ADP")

    u_c = reparametrisation.add_new_thermal_displacement_parameter(self.ind_ref)
    angle = reparametrisation.add(_.independent_scalar_parameter,
      value=self.angle_value*pi/180, variable=self.refine_angle)
    param = reparametrisation.add(
      _.shared_rotated_u_star,
      scatterer=scatterers[self.ind_atom],
      reference=u_c,
      direction=reparametrisation.find_direction(self.direction),
      angle=angle
    )
    reparametrisation.shared_Us[self.ind_atom] = u_c
    reparametrisation.asu_scatterer_parameters[self.ind_atom].u = param
    self.value = u_c
    self.angle = angle

class scalar_scaled_u(object):
  """ u_iso or u_star of a group of atoms all equal to their starting values
      times a scalar constant
  """

  def __init__(self, ind_sequence):
    self.indices = ind_sequence

  @property
  def constrained_parameters(self):
    return tuple((idx, "U") for idx in self.indices)

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    scalar = reparametrisation.add(_.independent_scalar_parameter,
      value=1., variable=True)
    for idx in self.indices:
      use_u_aniso = scatterers[idx].flags.use_u_aniso()
      if use_u_aniso:
        param = reparametrisation.add(
          _.scalar_scaled_u_star,
          scalar=scalar,
          scatterer=scatterers[idx])
      else:
        param = reparametrisation.add(
          _.scalar_scaled_u_iso,
          scalar=scalar,
          scatterer=scatterers[idx])
      # reparametrisation.shared_Us[idx] = param # not sure what to do with this
      reparametrisation.asu_scatterer_parameters[idx].u = param
    self.scalar = scalar

  def esd(self, ls):
    from math import sqrt
    cov_diag = ls.covariance_matrix().matrix_packed_u_diagonal()
    return sqrt(cov_diag[self.scalar.index])
