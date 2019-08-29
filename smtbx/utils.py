from __future__ import absolute_import, division, print_function
from cctbx import crystal, sgtbx
import libtbx

class connectivity_table(object):
  """ Bond connectivity tabulated """

  covalent_bond_tolerance = 0.5 # Angstrom
  conformer_indices=None
  sym_excl_indices=None
  radii=None
  # TODO: add possibility to fine tune connectivity by addition and deletion
  #       of individual bonds

  def __init__(self,
               structure,
               **kwds):
    from cctbx.eltbx import covalent_radii
    self.structure = structure
    libtbx.adopt_optional_init_args(self, kwds)
    max_r = 0
    for st in structure.scattering_type_registry().type_index_pairs_as_dict():
      r = 0
      if self.radii:
        r = self.radii.get(st, 0)
      if r == 0:
        r = covalent_radii.table(st).radius()
      if r > max_r: max_r = r
    self.structure = structure
    self.buffer_thickness = 2*max_r + self.covalent_bond_tolerance
    asu_mappings = structure.asu_mappings(
      buffer_thickness=self.buffer_thickness)
    self._pair_asu_table = crystal.pair_asu_table(asu_mappings)
    self._pair_asu_table_needs_updating = False
    if self.radii is None:
      self.radii = {}
    self._pair_asu_table.add_covalent_pairs(
      structure.scattering_types(),
      conformer_indices=self.conformer_indices,
      sym_excl_indices=self.sym_excl_indices,
      tolerance=self.covalent_bond_tolerance,
      radii=self.radii
    )
    self.pair_sym_table = self.pair_asu_table.extract_pair_sym_table()


  @property
  def pair_asu_table(self):
    if self._pair_asu_table_needs_updating:
      self._pair_asu_table = crystal.pair_asu_table(
        asu_mappings=self._pair_asu_table.asu_mappings())
      self._pair_asu_table.add_pair_sym_table(self.pair_sym_table)
      self.pair_sym_table = self._pair_asu_table.extract_pair_sym_table()
    return self._pair_asu_table

  def remove_bond(self, i_seq, j_seq, rt_mx_ji=sgtbx.rt_mx()):
    space_group = self.pair_asu_table.asu_mappings().space_group()
    r_den, t_den = space_group.r_den(), space_group.t_den()
    if j_seq < i_seq:
      i_seq, j_seq = j_seq, i_seq
      if not rt_mx_ji.is_unit_mx():
        rt_mx_ji = rt_mx_ji.inverse()
    if j_seq not in self.pair_sym_table[i_seq]:
      return
    for i, rt_mx in enumerate(self.pair_sym_table[i_seq][j_seq]):
      if (rt_mx.new_denominators(r_den, t_den) ==
          rt_mx_ji.new_denominators(r_den, t_den)):
        del self.pair_sym_table[i_seq][j_seq][i]
    self._pair_asu_table_needs_updating = True

  def add_bond(self, i_seq, j_seq, rt_mx_ji=sgtbx.rt_mx()):
    try:
      self.pair_asu_table.add_pair(i_seq, j_seq, rt_mx_ji)
    except RuntimeError:
      sites_frac = self.structure.sites_frac()
      sites = [sites_frac[i_seq], sites_frac[j_seq]]
      if not rt_mx_ji.is_unit_mx():
        sites[-1] = rt_mx_ji * sites[-1]
      d = self.structure.unit_cell().distance(sites[0], sites[1])
      self.pair_sym_table = self.pair_asu_table.extract_pair_sym_table()
      self._pair_asu_table = crystal.pair_asu_table(self.structure.asu_mappings(
        buffer_thickness=max(self.buffer_thickness, d)))
      self._pair_asu_table.add_pair_sym_table(self.pair_sym_table)
      self.pair_asu_table.add_pair(i_seq, j_seq, rt_mx_ji)
