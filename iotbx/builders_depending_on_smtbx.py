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

  def _checked_scatterer(self, index, what):
    """A scatterer by index, or an error saying which reference is dangling.

    Indexing the array directly gives "Index out of range" and nothing else -
    not the constraint, not the atom, not the file it came from. That happens
    whenever a constraint outlives the atom ordering it was built against: a
    RESI that takes a riding group's pivot into a residue and leaves the
    hydrogens outside is enough to do it, and the message named neither RESI
    nor an atom.
    """
    sc = self.structure.scatterers()
    if not (0 <= index < len(sc)):
      known = ", ".join(s.label for s in sc[:6])
      raise RuntimeError(
        "%s refers to scatterer %d, but the structure has %d "
        "(%s%s). A constraint or a restraint is left over from a different "
        "atom ordering - check any AFIX group split across a RESI, or an "
        "instruction naming an atom that has been renamed or deleted."
        % (what, index, len(sc), known, ", ..." if len(sc) > 6 else ""))
    return sc[index]

  def add_occupancy_pair_affine_constraint(self, scatterer_indices, linear_form):
    """ Add a constraint on the occupancies of a pair of scatterers that is
        affine, i.e. linear_form shall be ((a0, a1), b) and then
           a0*occ0 + a1*occ1 = b
        where (occ0, occ1) are the occupancies of the scatterers whose indices
        are given in `scatterer_indices`.
    """
    for i in scatterer_indices:
      self._checked_scatterer(i, "An occupancy constraint")
    self.constraints.append(
      constraints.occupancy.occupancy_pair_affine_constraint(scatterer_indices,
                                                             linear_form))

  def add_u_iso_proportional_to_pivot_u_eq(self,
                                           u_iso_scatterer_index,
                                           u_eq_scatterer_index,
                                           multiplier):
    sc_eq = self._checked_scatterer(
      u_eq_scatterer_index, "A riding-ADP constraint's pivot")
    sc_iso = self._checked_scatterer(
      u_iso_scatterer_index, "A riding-ADP constraint's dependent atom")
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
