from cctbx.crystal.find_best_cell import find_best_cell
from cctbx import uctbx
from cctbx import sgtbx

import sys

class symmetry(object):

  def __init__(self, unit_cell=None,
                     space_group_symbol=None,
                     space_group_info=None,
                     space_group=None,
                     assert_is_compatible_unit_cell=0001,
                     force_compatible_unit_cell=0001):
    assert [space_group_symbol, space_group_info, space_group].count(None) >= 2
    if (    unit_cell is not None
        and not isinstance(unit_cell, uctbx.ext.unit_cell)):
      unit_cell = uctbx.unit_cell(unit_cell)
    self._unit_cell = unit_cell
    self._space_group_info = space_group_info
    if (self._space_group_info is None):
      if (space_group_symbol is not None):
        self._space_group_info = sgtbx.space_group_info(space_group_symbol)
      elif (space_group is not None):
        self._space_group_info = sgtbx.space_group_info(group=space_group)

    if (self.unit_cell() is not None and self.space_group_info() is not None):
      if (assert_is_compatible_unit_cell):
        assert self.is_compatible_unit_cell(), \
          "Space group is incompatible with unit cell parameters."
      if (force_compatible_unit_cell):
        self._unit_cell = self.space_group().average_unit_cell(self._unit_cell)

  def _copy_constructor(self, other):
    self._unit_cell = other._unit_cell
    self._space_group_info = other._space_group_info

  def unit_cell(self):
    return self._unit_cell

  def space_group_info(self):
    return self._space_group_info

  def space_group(self):
    return self.space_group_info().group()

  def show_summary(self, f=sys.stdout):
    if (self.unit_cell() is None):
      print >> f, "Unit cell:", None
    else:
      self.unit_cell().show_parameters(f)
    if (self.space_group_info() is None):
      print >> f, "Space group:", None
    else:
      self.space_group_info().show_summary(f)

  def is_compatible_unit_cell(self):
    return self.space_group().is_compatible_unit_cell(self.unit_cell())

  def cell_equivalent_p1(self):
    return symmetry(self.unit_cell(), space_group_symbol="P 1")

  def change_basis(self, cb_op):
    if (isinstance(cb_op, str)):
      cb_op = sgtbx.change_of_basis_op(cb_op)
    return symmetry(
      unit_cell=cb_op.apply(self.unit_cell()),
      space_group_info=self.space_group_info().change_basis(cb_op))

  def primitive_setting(self):
    return self.change_basis(self.space_group().z2p_op())

  def change_of_basis_op_to_reference_setting(self):
    return self.space_group_info().type().cb_op()

  def as_reference_setting(self):
    return self.change_basis(self.change_of_basis_op_to_reference_setting())

  def change_of_basis_op_to_best_cell(self, angular_tolerance=None):
    return find_best_cell(self, angular_tolerance=angular_tolerance).cb_op()

  def best_cell(self, angular_tolerance=None):
    return self.change_basis(self.change_of_basis_op_to_best_cell(
      angular_tolerance=angular_tolerance))

  def change_of_basis_op_to_minimum_cell(self):
    z2p_op = self.space_group().z2p_op()
    r_inv = z2p_op.c_inv().r()
    p_cell = self.unit_cell().change_basis(r_inv.num(), r_inv.den())
    red = p_cell.minimum_reduction()
    p2n_op = sgtbx.change_of_basis_op(
      sgtbx.rt_mx(sgtbx.rot_mx(red.r_inv(), 1))).inverse()
    return p2n_op.new_denominators(z2p_op) * z2p_op

  def minimum_cell(self):
    return self.change_basis(self.change_of_basis_op_to_minimum_cell())

  def change_of_basis_op_to_niggli_cell(self):
    z2p_op = self.space_group().z2p_op()
    r_inv = z2p_op.c_inv().r()
    p_cell = self.unit_cell().change_basis(r_inv.num(), r_inv.den())
    red = p_cell.niggli_reduction()
    p2n_op = sgtbx.change_of_basis_op(
      sgtbx.rt_mx(sgtbx.rot_mx(red.r_inv().elems, 1))).inverse()
    return p2n_op.new_denominators(z2p_op) * z2p_op

  def niggli_cell(self):
    return self.change_basis(self.change_of_basis_op_to_niggli_cell())

  def patterson_symmetry(self):
    return symmetry(
      unit_cell=self.unit_cell(),
      space_group=self.space_group().build_derived_patterson_group())

  def is_patterson_symmetry(self):
    return self.space_group().build_derived_patterson_group() \
        == self.space_group()

  def join_symmetry(self, other_symmetry, force=00000):
    if (other_symmetry is None):
      return self
    if (force == 00000):
      strong = self
      weak = other_symmetry
    else:
      strong = other_symmetry
      weak = self
    unit_cell = strong.unit_cell()
    space_group_info = strong.space_group_info()
    if (unit_cell is None):
      unit_cell = weak.unit_cell()
    if (space_group_info is None):
      space_group_info = weak.space_group_info()
    return symmetry(
       unit_cell=unit_cell,
       space_group_info=space_group_info)

class special_position_settings(symmetry):

  def __init__(self, crystal_symmetry,
               min_distance_sym_equiv=0.5,
               u_star_tolerance=0,
               assert_is_positive_definite=00000,
               assert_min_distance_sym_equiv=0001):
    symmetry._copy_constructor(self, crystal_symmetry)
    self._min_distance_sym_equiv = min_distance_sym_equiv
    self._u_star_tolerance = u_star_tolerance
    self._assert_is_positive_definite = assert_is_positive_definite
    self._assert_min_distance_sym_equiv = assert_min_distance_sym_equiv

  def _copy_constructor(self, other):
    symmetry._copy_constructor(self, other)
    self._min_distance_sym_equiv = other._min_distance_sym_equiv
    self._u_star_tolerance = other._u_star_tolerance
    self._assert_is_positive_definite = other._assert_is_positive_definite
    self._assert_min_distance_sym_equiv = other._assert_min_distance_sym_equiv

  def min_distance_sym_equiv(self):
    return self._min_distance_sym_equiv

  def u_star_tolerance(self):
    return self._u_star_tolerance

  def assert_is_positive_definite(self):
    return self._assert_is_positive_definite

  def assert_min_distance_sym_equiv(self):
    return self._assert_min_distance_sym_equiv

  def site_symmetry(self, site):
    return sgtbx.site_symmetry(
      self.unit_cell(),
      self.space_group(),
      site,
      self.min_distance_sym_equiv(),
      self.assert_min_distance_sym_equiv())

  def sym_equiv_sites(self, site):
    return sgtbx.sym_equiv_sites(self.site_symmetry(site))

  def change_basis(self, cb_op):
    return special_position_settings(
      crystal_symmetry=symmetry.change_basis(self, cb_op),
      min_distance_sym_equiv=self.min_distance_sym_equiv(),
      u_star_tolerance=self.u_star_tolerance(),
      assert_is_positive_definite=self.assert_is_positive_definite(),
      assert_min_distance_sym_equiv=self.assert_min_distance_sym_equiv())
