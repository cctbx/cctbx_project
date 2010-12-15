from __future__ import division

from smtbx.refinement import constraints, least_squares
import smtbx.refinement.constraints.adp
import smtbx.refinement.constraints.geometrical.all

from iotbx.builders import \
     crystal_structure_builder, \
     restrained_crystal_structure_builder
import iotbx.constrained_parameters


class constrained_crystal_structure_builder(crystal_structure_builder):

  def __init__(self, *args, **kwds):
    super(constrained_crystal_structure_builder, self).__init__(*args, **kwds)
    self.constraints = []
    self.temperature_in_celsius = None

  def add_scatterer(self, scatterer, behaviour_of_variable, *args, **kwds):
    _ = iotbx.constrained_parameters
    crystal_structure_builder.add_scatterer(
      self, scatterer, behaviour_of_variable, *args, **kwds)
    if (scatterer.flags.use_u_iso()):
      b = behaviour_of_variable[4]
      if isinstance(b, tuple) and b[0] == _.constant_times_u_eq:
        self.constraints.append(
          constraints.adp.u_iso_proportional_to_pivot_u_eq(
            u_eq_scatterer_idx=b[2],
            u_iso_scatterer_idx=len(self.structure.scatterers()) - 1,
            multiplier=b[1]))

  def make_geometrical_constraint_type(self, constraint_name):
    return getattr(constraints.geometrical.all, constraint_name)

  def start_geometrical_constraint(self, type_,
                                   bond_length, rotating, stretching,
                                   pivot_relative_pos):
    self.first = len(self.structure.scatterers())

    self.current = type_(rotating=rotating,
                         stretching=stretching,
                         bond_length=bond_length,
                         pivot=self.first + pivot_relative_pos)

  def end_geometrical_constraint(self):
    last = len(self.structure.scatterers())
    self.current.constrained_site_indices = tuple(xrange(self.first, last))
    self.constraints.append(self.current)


class weighting_scheme_builder(object):

  def make_shelx_weighting_scheme(self, a, b, c=0, d=0, e=0, f=1/3):
    assert f == 1/3
    if c == 0 and d == 0 and e == 0:
      self.weighting_scheme = \
          least_squares.mainstream_shelx_weighting(a, b)
    else:
      self.weighting_scheme = \
          least_squares.shelx_weighting(a, b, c, d, e, f)

class constrained_restrained_crystal_structure_builder(
  constrained_crystal_structure_builder, restrained_crystal_structure_builder):
  pass


class weighted_constrained_restrained_crystal_structure_builder(
  weighting_scheme_builder, constrained_restrained_crystal_structure_builder):
  pass
