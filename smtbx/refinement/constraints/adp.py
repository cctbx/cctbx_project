import smtbx.refinement.constraints as _
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
    param = reparametrisation.add(
      _.u_iso_proportional_to_pivot_u_eq,
      pivot_u=reparametrisation.add_new_thermal_displacement_parameter(
        self.u_eq_scatterer_idx),
      scatterer = reparametrisation.structure.scatterers()[
        self.u_iso_scatterer_idx],
      multiplier=self.multiplier)
    reparametrisation.asu_scatterer_parameters[
      self.u_iso_scatterer_idx].u = param
