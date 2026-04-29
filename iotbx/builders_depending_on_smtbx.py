"""Tools for creating builders using smtbx
"""
from __future__ import absolute_import, division, print_function

from smtbx.refinement import constraints, least_squares
import smtbx.refinement.constraints.adp
import smtbx.refinement.constraints.geometrical.all
import smtbx.refinement.constraints.occupancy

from iotbx.builders import \
     crystal_structure_builder, \
     restrained_crystal_structure_builder
from six.moves import range


class constrained_crystal_structure_builder(crystal_structure_builder):

  def __init__(self, *args, **kwds):
    super(constrained_crystal_structure_builder, self).__init__(*args, **kwds)
    self.constraints = []
    self.temperature_in_celsius = None

  def add_occupancy_pair_affine_constraint(self, scatterer_indices, linear_form):
    """ Add a constraint on the occupancies of a pair of scatterers that is
        affine, i.e. linear_form shall be ((a0, a1), b) and then
           a0*occ0 + a1*occ1 = b
        where (occ0, occ1) are the occupancies of the scatterers whose indices
        are given in `scatterer_indices`.
    """
    self.constraints.append(
      constraints.occupancy.occupancy_pair_affine_constraint(scatterer_indices,
                                                             linear_form))

  def add_u_iso_proportional_to_pivot_u_eq(self,
                                           u_iso_scatterer_index,
                                           u_eq_scatterer_index,
                                           multiplier):
    sc = self.structure.scatterers()
    sc_eq = sc[u_eq_scatterer_index]
    sc_iso = sc[u_iso_scatterer_index]
    if sc_iso.flags.use_u_iso():
      self.constraints.append(
        constraints.adp.u_iso_proportional_to_pivot_u_eq(
          u_eq_scatterer_idx=u_eq_scatterer_index,
          u_iso_scatterer_idx=u_iso_scatterer_index,
          multiplier=multiplier))

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
    self.current.constrained_site_indices = tuple(range(self.first, last))
    self.constraints.append(self.current)


class weighting_scheme_builder(object):

  def make_shelx_weighting_scheme(self, a, b, c=0, d=0, e=0, f=1/3):
    assert f == 1/3, "Non-Wilsonian ShelX weighting not supported"
    if c == 0 and d == 0 and e == 0:
      self.weighting_scheme = \
          least_squares.mainstream_shelx_weighting(a, b)
    else:
      raise NotImplementedError(
        "ShelX weighting scheme with non-zero parameter c, d or e")

class weighted_constrained_restrained_crystal_structure_builder(
  constrained_crystal_structure_builder,
  restrained_crystal_structure_builder,
  weighting_scheme_builder):
  pass
