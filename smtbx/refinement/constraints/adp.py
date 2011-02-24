import smtbx.refinement.constraints as _
from smtbx.refinement.constraints import InvalidConstraint
import itertools

class u_iso_proportional_to_pivot_u_eq(object):
  """ u_iso of some scatterer constrained to be proportional to
      equivalent u_iso associated with adp of another scatterer
  """

  __slots__ = ('u_iso_scatterer_idx', 'u_eq_scatterer_idx', 'multiplier')

  def __init__(self, *args, **kwds):
    for attr, value in itertools.chain(
      itertools.izip(self.__slots__, args), kwds.iteritems()
      ):
      setattr(self, attr, value)

  def __eq__(self, other):
    try:
      for attr in self.__slots__:
        if getattr(self, attr) != getattr(other, attr): return False
      else:
        return True
    except AttributeError:
      return False

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
        scatterer = reparametrisation.structure.scatterers()[
          self.u_iso_scatterer_idx],
        multiplier=self.multiplier)
    reparametrisation.asu_scatterer_parameters[
      self.u_iso_scatterer_idx].u = param

class shared_u(object):
  """ u_iso or u_star of some scatterer constrained to be equal to
      u_iso or u_start of another scatterer
  """

  def __init__(self, ind_sequence):
    if len(ind_sequence) < 2:
      raise InvalidConstraint("at least two atoms are expected")
    self.indices = ind_sequence

  # any one scatterer can be used only in one of this constraint
  def __eq__(self, other):
    if self.indices != other.indices:  return False
    return True

  def add_to(self, reparametrisation):
    scatterers = reparametrisation.structure.scatterers()
    src_uses_u = scatterers[self.indices[0]].flags.use_u_aniso()
    for i in xrange(1, len(self.indices)):
      if scatterers[self.indices[i]].flags.use_u_aniso() != src_uses_u:
        raise InvalidConstraint(
          "mixing isotropic and anisotropic atoms is not allowed for shared ADP")

    u_c = reparametrisation.add_new_thermal_displacement_parameter(
      self.indices[0])
    for i in xrange(1, len(self.indices)):
      if src_uses_u:
        param = reparametrisation.add(
          _.shared_u_star,
          u_star=u_c,
          scatterer = reparametrisation.structure.scatterers()[
            self.indices[i]])
      else:
        param = reparametrisation.add(
          _.shared_u_iso,
          u_iso=u_c,
          scatterer = reparametrisation.structure.scatterers()[
            self.indices[i]])
      reparametrisation.shared_Us[self.indices[i]] = u_c
      reparametrisation.asu_scatterer_parameters[self.indices[i]].u = param
    self.value = u_c
