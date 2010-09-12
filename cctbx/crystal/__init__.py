from cctbx.array_family import flex

import boost.python
ext = boost.python.import_ext("cctbx_crystal_ext")
from cctbx_crystal_ext import *

from cctbx.crystal.find_best_cell import find_best_cell
from cctbx import sgtbx
from cctbx import uctbx

from scitbx.array_family import shared
from scitbx import stl
import scitbx.stl.set
import scitbx.stl.vector
import libtbx
from libtbx.utils import Keep
import sys

import scitbx.cubicle_neighbors
cubicles_max_memory_allocation_set(
  number_of_bytes=scitbx.cubicle_neighbors.cubicles_max_memory_allocation_get())

pair_sym_ops = sgtbx.stl_vector_rt_mx

pair_asu_j_sym_groups = scitbx.stl.vector.set_unsigned
pair_asu_j_sym_group = scitbx.stl.set.unsigned

class symmetry(object):

  def __init__(self, unit_cell=None,
                     space_group_symbol=None,
                     space_group_info=None,
                     space_group=None,
                     assert_is_compatible_unit_cell=True,
                     force_compatible_unit_cell=True):
    assert [space_group_symbol, space_group_info, space_group].count(None)>=2
    if (    unit_cell is not None
        and not isinstance(unit_cell, uctbx.ext.unit_cell)):
      unit_cell = uctbx.unit_cell(unit_cell)
    self._unit_cell = unit_cell
    self._space_group_info = space_group_info
    if (self._space_group_info is None):
      if (space_group_symbol is not None):
        self._space_group_info = sgtbx.space_group_info(
          symbol=space_group_symbol)
      elif (space_group is not None):
        if (isinstance(space_group, sgtbx.space_group)):
          self._space_group_info = sgtbx.space_group_info(group=space_group)
        else:
          self._space_group_info = sgtbx.space_group_info(symbol=space_group)
    if (self.unit_cell() is not None and self.space_group_info() is not None):
      if (assert_is_compatible_unit_cell):
        assert self.is_compatible_unit_cell(), \
          "Space group is incompatible with unit cell parameters."
      if (force_compatible_unit_cell):
        self._unit_cell = self.space_group().average_unit_cell(
          self._unit_cell)

  def _copy_constructor(self, other):
    self._unit_cell = other._unit_cell
    self._space_group_info = other._space_group_info

  def customized_copy(self, unit_cell=Keep, space_group_info=Keep):
    if (unit_cell is Keep): unit_cell = self._unit_cell
    if (space_group_info is Keep): space_group_info = self._space_group_info
    return symmetry(unit_cell=unit_cell, space_group_info=space_group_info)

  def unit_cell(self):
    return self._unit_cell

  def space_group_info(self):
    return self._space_group_info

  def space_group(self):
    sgi = self._space_group_info
    if (sgi is None): return None
    return sgi.group()

  def show_summary(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    if (self.unit_cell() is None):
      print >> f, prefix + "Unit cell:", None
    else:
      self.unit_cell().show_parameters(f=f, prefix=prefix+"Unit cell: ")
    if (self.space_group_info() is None):
      print >> f, prefix + "Space group:", None
    else:
      self.space_group_info().show_summary(f=f, prefix=prefix+"Space group: ")

  def is_similar_symmetry(self, other, relative_length_tolerance=0.01,
                                       absolute_angle_tolerance=1.):
    if (not self.unit_cell().is_similar_to(other.unit_cell(),
      relative_length_tolerance, absolute_angle_tolerance)): return False
    return self.space_group() == other.space_group()

  def is_compatible_unit_cell(self):
    return self.space_group().is_compatible_unit_cell(self.unit_cell())

  def cell_equivalent_p1(self):
    return symmetry(self.unit_cell(), space_group_symbol="P 1")

  def change_basis(self, cb_op):
    if (isinstance(cb_op, str)):
      cb_op = sgtbx.change_of_basis_op(cb_op)
    return symmetry(
      unit_cell=self.unit_cell().change_basis(cb_op),
      space_group_info=self.space_group_info().change_basis(cb_op))

  def change_of_basis_op_to_primitive_setting(self):
    return self.space_group().z2p_op()

  def primitive_setting(self):
    return self.change_basis(self.change_of_basis_op_to_primitive_setting())

  def change_of_basis_op_to_reference_setting(self):
    return self.space_group_info().type().cb_op()

  def as_reference_setting(self):
    return self.change_basis(self.change_of_basis_op_to_reference_setting())

  def change_of_basis_op_to_best_cell(self,
        angular_tolerance=None,
        best_monoclinic_beta=True):
    return find_best_cell(
      input_symmetry=self,
      angular_tolerance=angular_tolerance,
      best_monoclinic_beta=best_monoclinic_beta).cb_op()

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

  def change_of_basis_op_to_niggli_cell(self,
        relative_epsilon=None,
        iteration_limit=None):
    z2p_op = self.space_group().z2p_op()
    r_inv = z2p_op.c_inv().r()
    p_cell = self.unit_cell().change_basis(r_inv.num(), r_inv.den())
    red = p_cell.niggli_reduction(
      relative_epsilon=relative_epsilon,
      iteration_limit=iteration_limit)
    p2n_op = sgtbx.change_of_basis_op(
      sgtbx.rt_mx(sgtbx.rot_mx(red.r_inv().elems, 1))).inverse()
    return p2n_op.new_denominators(z2p_op) * z2p_op

  def niggli_cell(self,
        relative_epsilon=None,
        iteration_limit=None):
    return self.change_basis(self.change_of_basis_op_to_niggli_cell(
      relative_epsilon=relative_epsilon,
      iteration_limit=iteration_limit))

  def change_of_basis_op_to_inverse_hand(self):
    return self.space_group_info().type().change_of_hand_op()

  def inverse_hand(self):
    return self.change_basis(self.change_of_basis_op_to_inverse_hand())

  def reflection_intensity_symmetry(self, anomalous_flag):
    return symmetry(
      unit_cell=self.unit_cell(),
      space_group=self.space_group()
        .build_derived_reflection_intensity_group(
          anomalous_flag=anomalous_flag))

  def patterson_symmetry(self):
    return symmetry(
      unit_cell=self.unit_cell(),
      space_group=self.space_group().build_derived_patterson_group())

  def is_patterson_symmetry(self):
    return self.space_group().build_derived_patterson_group() \
        == self.space_group()

  def join_symmetry(self, other_symmetry, force=False):
    if (other_symmetry is None):
      return symmetry(
         unit_cell=self.unit_cell(),
         space_group_info=self.space_group_info())
    if (force == False):
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

  def subtract_continuous_allowed_origin_shifts(self, translation_cart):
    uc = self.unit_cell()
    return uc.orthogonalize(
      self.space_group_info().subtract_continuous_allowed_origin_shifts(
        translation_frac=uc.fractionalize(translation_cart)))

  def direct_space_asu(self):
    return self.space_group_info().direct_space_asu().define_metric(
      unit_cell=self.unit_cell())

  def gridding(self, d_min=None,
                     resolution_factor=None,
                     step=None,
                     symmetry_flags=None,
                     mandatory_factors=None,
                     max_prime=5,
                     assert_shannon_sampling=True):
    from cctbx import maptbx
    return maptbx.crystal_gridding(
      unit_cell=self.unit_cell(),
      d_min=d_min,
      resolution_factor=resolution_factor,
      step=step,
      symmetry_flags=symmetry_flags,
      space_group_info=self.space_group_info(),
      mandatory_factors=mandatory_factors,
      max_prime=max_prime,
      assert_shannon_sampling=assert_shannon_sampling)

  def asu_mappings(self, buffer_thickness, asu_is_inside_epsilon=None):
    import cctbx.crystal.direct_space_asu
    return direct_space_asu.asu_mappings(
      space_group=self.space_group(),
      asu=self.direct_space_asu().as_float_asu(
        is_inside_epsilon=asu_is_inside_epsilon),
      buffer_thickness=buffer_thickness)

  def average_u_cart(self, u_cart):
    from cctbx import adptbx
    return adptbx.u_star_as_u_cart(self.unit_cell(),
      self.space_group().average_u_star(
        adptbx.u_cart_as_u_star(self.unit_cell(), u_cart)))

  def average_b_cart(self, b_cart):
    return self.average_u_cart(u_cart=b_cart)

  def special_position_settings(self,
        min_distance_sym_equiv=0.5,
        u_star_tolerance=0,
        assert_min_distance_sym_equiv=True):
    return special_position_settings(
      crystal_symmetry=self,
      min_distance_sym_equiv=min_distance_sym_equiv,
      u_star_tolerance=u_star_tolerance,
      assert_min_distance_sym_equiv=assert_min_distance_sym_equiv)

  def build_miller_set(self, anomalous_flag, d_min, d_max=None):
    from cctbx import miller
    return miller.build_set(
      crystal_symmetry=self,
      anomalous_flag=anomalous_flag,
      d_min=d_min,
      d_max=d_max)

  def as_cif_block(self):
    import iotbx.cif
    return iotbx.cif.crystal_symmetry_as_cif_block(self).cif_block

def select_crystal_symmetry(
      from_command_line     = None,
      from_parameter_file   = None,
      from_coordinate_files = [None],
      from_reflection_files = [None]):
  """Select/construct a crystal symmetry from a list of various options"""
  tmp = [from_command_line, from_parameter_file]+from_coordinate_files \
        +from_reflection_files
  if tmp.count(None)==len(tmp):
    raise AssertionError("No unit cell and symmetry information supplied")

  result = symmetry(
    unit_cell=None,
    space_group_info=None)
  if (from_command_line is not None):
    result = result.join_symmetry(
      other_symmetry=from_command_line, force=False)
  if (from_parameter_file is not None):
    result = result.join_symmetry(
      other_symmetry=from_parameter_file, force=False)
  if (result.unit_cell() is None):
    for crystal_symmetry in from_reflection_files:
      if crystal_symmetry is not None:
        unit_cell = crystal_symmetry.unit_cell()
        if (unit_cell is not None):
          result = symmetry(
            unit_cell=unit_cell,
            space_group_info=result.space_group_info(),
            assert_is_compatible_unit_cell=False)
          break
  for crystal_symmetry in from_coordinate_files:
    if crystal_symmetry is not None:
      result = result.join_symmetry(
        other_symmetry=crystal_symmetry, force=False)
  if (result.space_group_info() is None):
    for crystal_symmetry in from_reflection_files:
      space_group_info = None
      if crystal_symmetry is not None:
        space_group_info = crystal_symmetry.space_group_info()
      if (space_group_info is not None):
        result = symmetry(
          unit_cell=result.unit_cell(),
          space_group_info=space_group_info,
          assert_is_compatible_unit_cell=False)
        break
  return result

def non_crystallographic_symmetry(
      sites_cart=None,
      sites_cart_min=None,
      sites_cart_max=None,
      buffer_layer=None,
      default_buffer_layer=0.5,
      min_unit_cell_length=0):
  return symmetry(
    unit_cell=uctbx.non_crystallographic_unit_cell(
      sites_cart=sites_cart,
      sites_cart_min=sites_cart_min,
      sites_cart_max=sites_cart_max,
      buffer_layer=buffer_layer,
      default_buffer_layer=default_buffer_layer,
      min_unit_cell_length=min_unit_cell_length),
    space_group=sgtbx.space_group())

class special_position_settings(symmetry):

  def __init__(self, crystal_symmetry,
               min_distance_sym_equiv=0.5,
               u_star_tolerance=0,
               assert_min_distance_sym_equiv=True):
    symmetry._copy_constructor(self, crystal_symmetry)
    self._min_distance_sym_equiv = min_distance_sym_equiv
    self._u_star_tolerance = u_star_tolerance
    self._assert_min_distance_sym_equiv = assert_min_distance_sym_equiv

  def _copy_constructor(self, other):
    symmetry._copy_constructor(self, other)
    self._min_distance_sym_equiv = other._min_distance_sym_equiv
    self._u_star_tolerance = other._u_star_tolerance
    self._assert_min_distance_sym_equiv = other._assert_min_distance_sym_equiv

  def min_distance_sym_equiv(self):
    return self._min_distance_sym_equiv

  def u_star_tolerance(self):
    return self._u_star_tolerance

  def assert_min_distance_sym_equiv(self):
    return self._assert_min_distance_sym_equiv

  def change_basis(self, cb_op):
    return special_position_settings(
      crystal_symmetry=symmetry.change_basis(self, cb_op),
      min_distance_sym_equiv=self.min_distance_sym_equiv(),
      u_star_tolerance=self.u_star_tolerance(),
      assert_min_distance_sym_equiv=self.assert_min_distance_sym_equiv())

  def site_symmetry(self, site=None, site_cart=None):
    assert [site, site_cart].count(None) == 1
    if (site_cart is not None):
      site = self.unit_cell().fractionalize(site_cart)
    return sgtbx.site_symmetry(
      self.unit_cell(),
      self.space_group(),
      site,
      self.min_distance_sym_equiv(),
      self.assert_min_distance_sym_equiv())

  def sym_equiv_sites(self, site):
    return sgtbx.sym_equiv_sites(self.site_symmetry(site))

  def site_symmetry_table(self,
        sites_frac=None,
        sites_cart=None,
        unconditional_general_position_flags=None):
    assert (sites_frac is None) != (sites_cart is None)
    if (sites_frac is None):
      sites_frac = self.unit_cell().fractionalize(sites_cart=sites_cart)
    result = sgtbx.site_symmetry_table()
    result.process(
      unit_cell=self.unit_cell(),
      space_group=self.space_group(),
      original_sites_frac=sites_frac,
      unconditional_general_position_flags=
        unconditional_general_position_flags,
      min_distance_sym_equiv=self.min_distance_sym_equiv(),
      assert_min_distance_sym_equiv=self.assert_min_distance_sym_equiv())
    return result

  def asu_mappings(self,
        buffer_thickness,
        sites_frac=None,
        sites_cart=None,
        site_symmetry_table=None,
        asu_is_inside_epsilon=None):
    asu_mappings = symmetry.asu_mappings(self,
      buffer_thickness=buffer_thickness,
      asu_is_inside_epsilon=asu_is_inside_epsilon)
    if (sites_frac is not None or sites_cart is not None):
      assert sites_frac is None or sites_cart is None
      if (sites_frac is None):
        sites_frac = self.unit_cell().fractionalize(sites_cart=sites_cart)
      if (site_symmetry_table is None):
        site_symmetry_table = self.site_symmetry_table(sites_frac=sites_frac)
      asu_mappings.process_sites_frac(
        original_sites=sites_frac,
        site_symmetry_table=site_symmetry_table)
    return asu_mappings

  def pair_generator(self,
        distance_cutoff,
        sites_frac=None,
        sites_cart=None,
        site_symmetry_table=None,
        asu_mappings_buffer_thickness=None,
        asu_is_inside_epsilon=None,
        minimal=False):
    assert sites_frac is not None or sites_cart is not None
    if (asu_mappings_buffer_thickness is None):
        asu_mappings_buffer_thickness = distance_cutoff
    asu_mappings = self.asu_mappings(
      buffer_thickness=asu_mappings_buffer_thickness,
      sites_frac=sites_frac,
      sites_cart=sites_cart,
      site_symmetry_table=site_symmetry_table,
      asu_is_inside_epsilon=asu_is_inside_epsilon)
    return neighbors_fast_pair_generator(
      asu_mappings=asu_mappings,
      distance_cutoff=distance_cutoff,
      minimal=minimal)

  def pair_asu_table(self,
        distance_cutoff,
        sites_frac=None,
        sites_cart=None,
        site_symmetry_table=None,
        asu_mappings_buffer_thickness=None,
        asu_is_inside_epsilon=None,
        min_cubicle_edge=5,
        distance_cutoff_epsilon=None):
    assert sites_frac is not None or sites_cart is not None
    if (asu_mappings_buffer_thickness is None):
        asu_mappings_buffer_thickness = distance_cutoff
    asu_mappings = self.asu_mappings(
      buffer_thickness=asu_mappings_buffer_thickness,
      sites_frac=sites_frac,
      sites_cart=sites_cart,
      site_symmetry_table=site_symmetry_table,
      asu_is_inside_epsilon=asu_is_inside_epsilon)
    result = pair_asu_table(asu_mappings=asu_mappings)
    if (distance_cutoff_epsilon is None):
        distance_cutoff_epsilon = asu_mappings.asu().is_inside_epsilon()
    result.add_all_pairs(
      distance_cutoff=distance_cutoff,
      min_cubicle_edge=min_cubicle_edge,
      epsilon=distance_cutoff_epsilon)
    return result

  def incremental_pairs(self,
        distance_cutoff,
        asu_is_inside_epsilon=None,
        asu_mappings_buffer_thickness=-1,
        cubicle_epsilon=-1):
    result = incremental_pairs(
      space_group=self.space_group(),
      asu=self.direct_space_asu().as_float_asu(
        is_inside_epsilon=asu_is_inside_epsilon),
      distance_cutoff=distance_cutoff,
      asu_mappings_buffer_thickness=asu_mappings_buffer_thickness,
      cubicle_epsilon=cubicle_epsilon)
    result.min_distance_sym_equiv = self._min_distance_sym_equiv
    result.assert_min_distance_sym_equiv = self._assert_min_distance_sym_equiv
    return result

  def site_cluster_analysis(self,
        min_distance=None,
        min_cross_distance=None,
        min_self_distance=None,
        general_positions_only=False,
        estimated_reduction_factor=4,
        asu_is_inside_epsilon=None,
        asu_mappings_buffer_thickness=-1,
        min_cubicle_edge=5,
        cubicle_epsilon=-1):
    if (min_cross_distance is None): min_cross_distance = min_distance
    if (min_self_distance is None): min_self_distance = min_distance
    if (min_self_distance is None): min_self_distance = min_cross_distance
    assert min_cross_distance is not None
    assert min_self_distance is not None
    result = site_cluster_analysis(
      space_group=self.space_group(),
      asu=self.direct_space_asu().as_float_asu(
        is_inside_epsilon=asu_is_inside_epsilon),
      min_cross_distance=min_cross_distance,
      min_self_distance=min_self_distance,
      general_positions_only=general_positions_only,
      estimated_reduction_factor=estimated_reduction_factor,
      asu_mappings_buffer_thickness=asu_mappings_buffer_thickness,
      min_cubicle_edge=min_cubicle_edge,
      cubicle_epsilon=cubicle_epsilon)
    result.min_distance_sym_equiv = self._min_distance_sym_equiv
    result.assert_min_distance_sym_equiv = self._assert_min_distance_sym_equiv
    return result

def correct_special_position(
      crystal_symmetry,
      special_op,
      site_frac=None,
      site_cart=None,
      site_label=None,
      tolerance=1,
      error_message="Corrupt gradient calculations."):
  """
  During refinement it is essential to reset special positions
  because otherwise rounding error accumulate over many cycles.
  """
  assert (site_frac is None) != (site_cart is None)
  unit_cell = crystal_symmetry.unit_cell()
  if (site_frac is None):
    site_frac = unit_cell.fractionalize(site_cart)
  site_special_frac = special_op * site_frac
  distance_moved = unit_cell.distance(site_special_frac, site_frac)
  if (distance_moved > tolerance):
    error_message += "\n  unit_cell: %s" % str(unit_cell)
    error_message += "\n  space_group_info: %s" % str(crystal_symmetry.space_group_info())
    error_message += "\n  special_op: %s" % str(special_op)
    if (site_label is not None):
      error_message += "\n  site_label: %s" % site_label
    error_message += "\n  site_frac: %s" % str(site_frac)
    error_message += "\n  site_special_frac: %s" % str(site_special_frac)
    error_message += "\n  distance_moved: %g" % distance_moved
    error_message += "\n  ****** This is a very critical error. ******"
    error_message += "\n  PLEASE send this output to"
    error_message += "\n"
    error_message += "\n    cctbx@cci.lbl.gov"
    error_message += "\n"
    error_message += "\n  to help us resolve the problem."
    error_message += "\n  Thank you in advance!"
    raise AssertionError(error_message)
  if (site_cart is None):
    return site_special_frac
  return unit_cell.orthogonalize(site_special_frac)

class _pair_asu_table(boost.python.injector, pair_asu_table):

  def as_nested_lists(self):
    result = []
    for i_seq, j_seq_dict in enumerate(self.table()):
      i_seq_list = [i_seq]
      for j_seq,j_sym_group in j_seq_dict.items():
        j_seq_list = [j_seq]
        for j_syms in j_sym_group:
          j_seq_list.append(list(j_syms))
        i_seq_list.append(j_seq_list)
      result.append(i_seq_list)
    return result

  def show(self, f=None, site_labels=None):
    if (f is None): f = sys.stdout
    if (site_labels is None):
      for i_seq, j_seq_dict in enumerate(self.table()):
        print >> f, "i_seq:", i_seq
        for j_seq,j_sym_group in j_seq_dict.items():
          print >> f, "  j_seq:", j_seq
          for j_syms in j_sym_group:
            print >> f, "    j_syms:", list(j_syms)
    else:
      assert len(site_labels) == self.table().size()
      for i_seq, j_seq_dict in enumerate(self.table()):
        print >> f, "%s(%d)" % (site_labels[i_seq], i_seq)
        for j_seq,j_sym_group in j_seq_dict.items():
          print >> f, "  %s(%d)" % (site_labels[j_seq], j_seq)
          for j_syms in j_sym_group:
            print >> f, "    j_syms:", list(j_syms)

  def show_distances(self,
        site_labels=None,
        sites_frac=None,
        sites_cart=None,
        show_cartesian=False,
        keep_pair_asu_table=False,
        out=None):
    return show_distances(
      pair_asu_table=self,
      site_labels=site_labels,
      sites_frac=sites_frac,
      sites_cart=sites_cart,
      show_cartesian=show_cartesian,
      keep_pair_asu_table=keep_pair_asu_table,
      out=out)

  def show_angles(self,
        site_labels=None,
        sites_frac=None,
        sites_cart=None,
        keep_pair_asu_table=False,
        out=None):
    return show_angles(
      pair_asu_table=self,
      site_labels=site_labels,
      sites_frac=sites_frac,
      sites_cart=sites_cart,
      keep_pair_asu_table=keep_pair_asu_table,
      out=out)

class show_distances(object):

  def __init__(self,
        pair_asu_table,
        site_labels=None,
        sites_frac=None,
        sites_cart=None,
        show_cartesian=False,
        keep_pair_asu_table=False,
        out=None):
    assert [sites_frac, sites_cart].count(None) == 1
    if (out is None): out = sys.stdout
    if (keep_pair_asu_table):
      self.pair_asu_table = pair_asu_table
    else:
      self.pair_asu_table = None
    self.distances = flex.double()
    self.pair_counts = flex.size_t()
    asu_mappings = pair_asu_table.asu_mappings()
    unit_cell = asu_mappings.unit_cell()
    if (sites_frac is None):
      sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
    if (site_labels is None):
      label_len = len("%d" % (sites_frac.size()+1))
      label_fmt = "site_%%0%dd" % label_len
      label_len += 5
    else:
      label_len = 1
      for label in site_labels:
        label_len = max(label_len, len(label))
      label_fmt = "%%-%ds" % (label_len+1)
    for i_seq,asu_dict in enumerate(pair_asu_table.table()):
      rt_mx_i_inv = asu_mappings.get_rt_mx(i_seq, 0).inverse()
      site_frac_i = sites_frac[i_seq]
      pair_count = 0
      dists = flex.double()
      j_seq_i_group = []
      for j_seq,j_sym_groups in asu_dict.items():
        site_frac_j = sites_frac[j_seq]
        for i_group,j_sym_group in enumerate(j_sym_groups):
          pair_count += j_sym_group.size()
          j_sym = j_sym_group[0]
          rt_mx_ji = rt_mx_i_inv.multiply(asu_mappings.get_rt_mx(j_seq, j_sym))
          distance = unit_cell.distance(site_frac_i, rt_mx_ji * site_frac_j)
          dists.append(distance)
          j_seq_i_group.append((j_seq,i_group))
      if (site_labels is None):
        s = label_fmt % (i_seq+1)
      else:
        s = label_fmt % site_labels[i_seq]
      s += " pair count: %3d" % pair_count
      if (show_cartesian):
        formatted_site = [" %7.2f" % x
          for x in unit_cell.orthogonalize(site_frac_i)]
      else:
        formatted_site = [" %7.4f" % x for x in site_frac_i]
      print >> out, ("%%-%ds" % (label_len+23)) % s, \
        "<<"+",".join(formatted_site)+">>"
      permutation = flex.sort_permutation(data=dists)
      for j_seq,i_group in flex.select(j_seq_i_group, permutation):
        site_frac_j = sites_frac[j_seq]
        j_sym_groups = asu_dict[j_seq]
        j_sym_group = j_sym_groups[i_group]
        for i_j_sym,j_sym in enumerate(j_sym_group):
          rt_mx_ji = rt_mx_i_inv.multiply(
            asu_mappings.get_rt_mx(j_seq, j_sym))
          site_frac_ji = rt_mx_ji * site_frac_j
          distance = unit_cell.distance(site_frac_i, site_frac_ji)
          self.distances.append(distance)
          if (site_labels is None):
            print >> out, " ", label_fmt % (j_seq+1) + ":",
          else:
            print >> out, " ", label_fmt % (site_labels[j_seq] + ":"),
          print >> out, "%8.4f" % distance,
          if (i_j_sym != 0):
            s = "sym. equiv."
          else:
            s = "           "
          if (show_cartesian):
            formatted_site = [" %7.2f" % x
              for x in unit_cell.orthogonalize(site_frac_ji)]
          else:
            formatted_site = [" %7.4f" % x for x in site_frac_ji]
          s += " (" + ",".join(formatted_site) +")"
          print >> out, s
      if (pair_count == 0):
        print >> out, "  no neighbors"
      self.pair_counts.append(pair_count)

class show_angles(object):

  def __init__(self,
        pair_asu_table,
        site_labels=None,
        sites_frac=None,
        sites_cart=None,
        show_cartesian=False,
        keep_pair_asu_table=False,
        out=None):
    assert [sites_frac, sites_cart].count(None) == 1
    if (out is None): out = sys.stdout
    if (keep_pair_asu_table):
      self.pair_asu_table = pair_asu_table
    else:
      self.pair_asu_table = None
    self.distances = flex.double()
    self.angles = flex.double()
    rt_mxs = []
    self.pair_counts = flex.size_t()
    asu_mappings = pair_asu_table.asu_mappings()
    unit_cell = asu_mappings.unit_cell()
    if (sites_frac is None):
      sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
    if (site_labels is None):
      label_len = len("%d" % (sites_frac.size()+1))
      label_fmt = "site_%%0%dd" % label_len
      label_len += 5
    else:
      label_len = 1
      for label in site_labels:
        label_len = max(label_len, len(label))
      label_fmt = "%%-%ds" % (label_len+4)
      label_fmt *= 3
    ## angle is formed by j_seq-i_seq-k_seq
    for i_seq,asu_dict in enumerate(pair_asu_table.table()):
      rt_mx_i_inv = asu_mappings.get_rt_mx(i_seq, 0).inverse()
      site_frac_i = sites_frac[i_seq]
      angles = flex.double()
      for j_seq,j_sym_groups in asu_dict.items():
        site_frac_j = sites_frac[j_seq]
        for j_sym_group in j_sym_groups:
          for i_j_sym,j_sym in enumerate(j_sym_group):
            rt_mx_ji = rt_mx_i_inv.multiply(
              asu_mappings.get_rt_mx(j_seq, j_sym))
            site_frac_ji = rt_mx_ji * site_frac_j
            for k_seq, k_sym_groups in asu_dict.items():
              if k_seq == j_seq and j_sym_group.size() <= 1: continue
              if k_seq > j_seq: continue
              site_frac_k = sites_frac[k_seq]
              for k_sym_group in k_sym_groups:
                for i_k_sym,k_sym in enumerate(k_sym_group):
                  if j_seq == k_seq and i_j_sym <= i_k_sym: continue
                  if i_seq == k_seq and i_k_sym == 0: continue
                  rt_mx_ki = rt_mx_i_inv.multiply(
                    asu_mappings.get_rt_mx(k_seq, k_sym))
                  site_frac_ki = rt_mx_ki * site_frac_k
                  angle = unit_cell.angle(site_frac_ji, site_frac_i, site_frac_ki)
                  if angle is None: continue
                  self.angles.append(angle)
                  if (site_labels is None):
                    s = label_fmt % (j_seq+1) + ":"
                  else:
                    i_label = site_labels[i_seq]
                    j_label = site_labels[j_seq]
                    k_label = site_labels[k_seq]
                    if i_j_sym != 0:
                      if rt_mx_ji in rt_mxs:
                        j = rt_mxs.index(rt_mx_ji) + 1
                      else:
                        rt_mxs.append(rt_mx_ji)
                        j = len(rt_mxs)
                      j_label += "*%s" %j
                    if i_k_sym != 0:
                      if rt_mx_ki in rt_mxs:
                        k = rt_mxs.index(rt_mx_ki) + 1
                      else:
                        rt_mxs.append(rt_mx_ki)
                        k = len(rt_mxs)
                      k_label += "*%s" %k
                    s = label_fmt % (j_label, i_label, k_label)
                  s += " %6.2f" % angle
                  print >> out, s
    for i, rt_mx in enumerate(rt_mxs):
      print >> out, "*%s" %(i+1),
      print >> out, rt_mx

class sym_pair(libtbx.slots_getstate_setstate):

  __slots__ = ["i_seq", "j_seq", "rt_mx_ji"]

  def __init__(self, i_seq, j_seq, rt_mx_ji):
    self.i_seq = i_seq
    self.j_seq = j_seq
    self.rt_mx_ji = rt_mx_ji

  def i_seqs(self):
    return (self.i_seq, self.j_seq)

class _pair_sym_table(boost.python.injector, pair_sym_table):

  def iterator(self):
    for i_seq,pair_sym_dict in enumerate(self):
      for j_seq,sym_ops in pair_sym_dict.items():
        for rt_mx_ji in sym_ops:
          yield sym_pair(i_seq=i_seq, j_seq=j_seq, rt_mx_ji=rt_mx_ji)

  def show(self, f=None, site_labels=None):
    if (f is None): f = sys.stdout
    if (site_labels is None):
      for i_seq,pair_sym_dict in enumerate(self):
        print >> f, "i_seq:", i_seq
        for j_seq,sym_ops in pair_sym_dict.items():
          print >> f, "  j_seq:", j_seq
          for sym_op in sym_ops:
            print >> f, "   ", sym_op
    else:
      for i_seq,pair_sym_dict in enumerate(self):
        print >> f, "%s(%d)" % (site_labels[i_seq], i_seq)
        for j_seq,sym_ops in pair_sym_dict.items():
          print >> f, "  %s(%d)" % (site_labels[j_seq], j_seq)
          for sym_op in sym_ops:
            print >> f, "   ", sym_op

  def number_of_pairs_involving_symmetry(self):
    result = 0
    for i_seq,pair_sym_dict in enumerate(self):
      for j_seq,sym_ops in pair_sym_dict.items():
        for sym_op in sym_ops:
          if (not sym_op.is_unit_mx()):
            result += 1
    return result

  def simple_edge_list(self):
    result = []
    for i_seq,pair_sym_dict in enumerate(self):
      for j_seq,sym_ops in pair_sym_dict.items():
        for sym_op in sym_ops:
          if (sym_op.is_unit_mx()):
            result.append((i_seq, j_seq))
            break
    return result

  def full_simple_connectivity(self):
    result = shared.stl_set_unsigned(self.size())
    for i_seq,pair_sym_dict in enumerate(self):
      for j_seq,sym_ops in pair_sym_dict.items():
        for sym_op in sym_ops:
          if (sym_op.is_unit_mx()):
            result[i_seq].insert(j_seq)
            result[j_seq].insert(i_seq)
            break
    return result

  def is_paired(self, i_seq):
    if (self[i_seq].size() != 0): return True
    for pair_sym_dict in self:
      if (i_seq in pair_sym_dict): return True
    return False

class _clustering_mix_in(object):

  def sites_cart(self):
    return self.special_position_settings.unit_cell().orthogonalize(
      sites_frac=self.sites_frac)

  def tidy_index_groups_in_place(self):
    cluster_sizes = flex.size_t()
    for ig in self.index_groups:
      cluster_sizes.append(len(ig))
    permutation = flex.sort_permutation(cluster_sizes, reverse=True)
    sorted_index_groups = flex.select(
      sequence=self.index_groups,
      permutation=permutation)
    index_groups = []
    for i,ig in enumerate(sorted_index_groups):
      if (len(ig) == 0): break
      ig = flex.size_t(ig)
      index_groups.append(ig.select(flex.sort_permutation(data=ig)))
    self.index_groups = index_groups
    return self

def _priorities_based_on_cluster_size(scores, index_groups):
  cluster_sizes = flex.size_t()
  for ig in index_groups:
    cluster_sizes.append(len(ig))
  permutation = flex.sort_permutation(cluster_sizes, reverse=True)
  sorted_index_groups = flex.select(
    sequence=index_groups,
    permutation=permutation)
  priorities = flex.size_t()
  for i,ig in enumerate(sorted_index_groups):
    ig = flex.size_t(ig)
    priorities.extend(ig.select(
     flex.sort_permutation(data=scores.select(ig), reverse=True)))
  return priorities

class incremental_clustering(_clustering_mix_in):

  def __init__(self,
        special_position_settings,
        sites_cart,
        distance_cutoffs,
        scores=None,
        discard_special_positions=False,
        discard_not_strictly_inside_asu=False,
        initial_required_cluster_size=0):
    self.special_position_settings = special_position_settings
    self.distance_cutoffs = distance_cutoffs
    unit_cell = special_position_settings.unit_cell()
    sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
    site_symmetry_table = special_position_settings.site_symmetry_table(
      sites_frac=sites_frac)
    if (not discard_special_positions):
      for i_seq in site_symmetry_table.special_position_indices():
        sites_frac[i_seq] = site_symmetry_table.get(i_seq).special_op() \
                          * sites_frac[i_seq]
    else:
      selection = ~flex.bool(
        sites_frac.size(),
        site_symmetry_table.special_position_indices())
      site_symmetry_table = site_symmetry_table.select(selection)
      assert site_symmetry_table.n_special_positions() == 0
      sites_frac = sites_frac.select(selection)
      if (scores is not None):
        scores = scores.select(selection)
    if (discard_not_strictly_inside_asu):
      selection = special_position_settings \
        .direct_space_asu() \
        .as_float_asu().is_inside_frac(sites_frac=sites_frac)
      site_symmetry_table = site_symmetry_table.select(selection)
      sites_frac = sites_frac.select(selection)
      if (scores is not None):
        scores = scores.select(selection)
    n_sites = sites_frac.size()
    assignments = flex.size_t(xrange(n_sites))
    n_clusters = n_sites
    index_groups = [[i_seq] for i_seq in xrange(n_sites)]
    symmetry_ops = [special_position_settings.space_group()(0)] * n_sites
    get_distance = unit_cell.distance
    for i_distance_cutoff,distance_cutoff in enumerate(distance_cutoffs):
      if (n_clusters == 1): break
      asu_mappings = special_position_settings.asu_mappings(
        buffer_thickness=distance_cutoff,
        sites_frac=sites_frac,
        site_symmetry_table=site_symmetry_table)
      pair_asu_tab = pair_asu_table(asu_mappings=asu_mappings)
      pair_asu_tab.add_all_pairs(distance_cutoff=distance_cutoff)
      if (scores is None):
        scores = pair_asu_tab.pair_counts()
      if (n_clusters == n_sites):
        priorities = flex.sort_permutation(data=scores, reverse=True)
      else:
        priorities = _priorities_based_on_cluster_size(
          scores=scores,
          index_groups=index_groups)
      for i_seq in priorities:
        i_cluster = assignments[i_seq]
        if (    i_distance_cutoff != 0
            and i_distance_cutoff != len(distance_cutoffs)-1
            and len(index_groups[i_cluster])<initial_required_cluster_size):
          continue
        rt_mx_ip = asu_mappings.get_rt_mx(i_seq, 0)
        site_ip = rt_mx_ip * sites_frac[i_seq]
        j_seq_dict = pair_asu_tab.table()[i_seq]
        for j_seq,j_sym_group in j_seq_dict.items():
          j_cluster = assignments[j_seq]
          if (j_cluster == i_cluster): continue
          if (    i_distance_cutoff != 0
              and i_distance_cutoff != len(distance_cutoffs)-1
              and len(index_groups[j_cluster])<initial_required_cluster_size):
            continue
          smallest_distance = distance_cutoff*(1+1.e-6)
          best_rt_mx_jp = None
          for j_syms in j_sym_group:
            j_sym = j_syms[0]
            rt_mx_jp = asu_mappings.get_rt_mx(j_seq, j_sym)
            site_jp = rt_mx_jp * sites_frac[j_seq]
            distance = get_distance(site_ip, site_jp)
            if (smallest_distance > distance):
              smallest_distance = distance
              best_rt_mx_jp = rt_mx_jp
          assert best_rt_mx_jp is not None
          rt_mx_jp = best_rt_mx_jp
          if (   len(index_groups[i_cluster])
              >= len(index_groups[j_cluster])):
            rt_mx_c = symmetry_ops[i_seq].multiply(
              rt_mx_ip.inverse().multiply(
                rt_mx_jp.multiply(
                  symmetry_ops[j_seq].inverse())))
            for c_seq in index_groups[j_cluster]:
              index_groups[i_cluster].append(c_seq)
              assignments[c_seq] = i_cluster
              symmetry_ops[c_seq] = rt_mx_c.multiply(symmetry_ops[c_seq])
            index_groups[j_cluster] = []
          else:
            rt_mx_c = symmetry_ops[j_seq].multiply(
              rt_mx_jp.inverse().multiply(
                rt_mx_ip.multiply(
                  symmetry_ops[i_seq].inverse())))
            for c_seq in index_groups[i_cluster]:
              index_groups[j_cluster].append(c_seq)
              assignments[c_seq] = j_cluster
              symmetry_ops[c_seq] = rt_mx_c.multiply(symmetry_ops[c_seq])
            index_groups[i_cluster] = []
            i_cluster = j_cluster
          n_clusters -= 1
    assert n_clusters > 0
    for i_seq,site_frac in enumerate(sites_frac):
      sites_frac[i_seq] = symmetry_ops[i_seq] * site_frac
    self.sites_frac = sites_frac
    self.index_groups = index_groups

class distance_based_clustering(_clustering_mix_in):

  def __init__(self,
        special_position_settings,
        sites_cart,
        distance_cutoff,
        discard_special_positions=False):
    self.special_position_settings = special_position_settings
    self.distance_cutoff = distance_cutoff
    unit_cell = special_position_settings.unit_cell()
    sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
    site_symmetry_table = special_position_settings.site_symmetry_table(
      sites_frac=sites_frac)
    if (not discard_special_positions):
      for i_seq in site_symmetry_table.special_position_indices():
        sites_frac[i_seq] = site_symmetry_table.get(i_seq).special_op() \
                          * sites_frac[i_seq]
    else:
      selection = ~flex.bool(
        sites_frac.size(),
        site_symmetry_table.special_position_indices())
      site_symmetry_table = site_symmetry_table.select(selection)
      assert site_symmetry_table.n_special_positions() == 0
      sites_frac = sites_frac.select(selection)
    n_sites = sites_frac.size()
    assignments = flex.size_t(xrange(n_sites))
    n_clusters = n_sites
    index_groups = [[i_seq] for i_seq in xrange(n_sites)]
    symmetry_ops = [special_position_settings.space_group()(0)] * n_sites
    asu_mappings = special_position_settings.asu_mappings(
      buffer_thickness=distance_cutoff,
      sites_frac=sites_frac,
      site_symmetry_table=site_symmetry_table)
    pair_asu_tab = pair_asu_table(asu_mappings=asu_mappings)
    pair_asu_tab.add_all_pairs(distance_cutoff=distance_cutoff)
    pair_sym_tab = pair_asu_tab.extract_pair_sym_table()
    sym_pairs = []
    distances = flex.double()
    get_distance = unit_cell.distance
    for pair in pair_sym_tab.iterator():
      site_i = sites_frac[pair.i_seq]
      site_j = sites_frac[pair.j_seq]
      site_ji = pair.rt_mx_ji * site_j
      distance = get_distance(site_i, site_ji)
      assert distance <= distance_cutoff * (1+1.e-6), distance
      sym_pairs.append(pair)
      distances.append(distance)
    priorities = flex.sort_permutation(data=distances)
    for i_pair in priorities:
      if (n_clusters == 1): break
      pair = sym_pairs[i_pair]
      i_seq = pair.i_seq
      j_seq = pair.j_seq
      i_cluster = assignments[i_seq]
      j_cluster = assignments[j_seq]
      if (i_cluster == j_cluster): continue
      if (   len(index_groups[i_cluster])
          >= len(index_groups[j_cluster])):
        rt_mx_c = symmetry_ops[i_seq].multiply(
          pair.rt_mx_ji.multiply(
            symmetry_ops[j_seq].inverse()))
        for c_seq in index_groups[j_cluster]:
          index_groups[i_cluster].append(c_seq)
          assignments[c_seq] = i_cluster
          symmetry_ops[c_seq] = rt_mx_c.multiply(symmetry_ops[c_seq])
        index_groups[j_cluster] = []
      else:
        rt_mx_c = symmetry_ops[j_seq].multiply(
          pair.rt_mx_ji.inverse().multiply(
            symmetry_ops[i_seq].inverse()))
        for c_seq in index_groups[i_cluster]:
          index_groups[j_cluster].append(c_seq)
          assignments[c_seq] = j_cluster
          symmetry_ops[c_seq] = rt_mx_c.multiply(symmetry_ops[c_seq])
        index_groups[i_cluster] = []
      n_clusters -= 1
      assert abs(get_distance(
        symmetry_ops[i_seq]*sites_frac[i_seq],
        symmetry_ops[j_seq]*sites_frac[j_seq]) - distances[i_pair]) < 1.e-6
    assert n_clusters > 0
    for i_seq,site_frac in enumerate(sites_frac):
      sites_frac[i_seq] = symmetry_ops[i_seq] * site_frac
    self.sites_frac = sites_frac
    self.index_groups = index_groups

def cluster_erosion(sites_cart, box_size, fraction_to_be_retained):
  from scitbx import matrix
  from libtbx.math_utils import iceil
  lower_left = matrix.col(sites_cart.min())
  upper_right = matrix.col(sites_cart.max())
  sites_span = upper_right - lower_left
  n_boxes = matrix.col([iceil(x/box_size) for x in sites_span])
  box_span = n_boxes * box_size
  sites_center = (lower_left + upper_right) / 2
  box_origin = sites_center - box_span / 2
  box_coordinates = (sites_cart - box_origin) * (1./box_size)
  box_grid = flex.grid(n_boxes)
  box_counts = flex.size_t(box_grid.size_1d(), 0)
  box_members = [flex.size_t() for i in xrange(box_counts.size())]
  for i_seq,x in enumerate(box_coordinates):
    box_index_3d = [max(0, min(int(i), n)) for i,n in zip(x, n_boxes)]
    box_index_1d = box_grid(box_index_3d)
    box_counts[box_index_1d] += 1
    box_members[box_index_1d].append(i_seq)
  perm = flex.sort_permutation(data=box_counts, reverse=True)
  box_counts_sorted = box_counts.select(perm)
  threshold = sites_cart.size() * fraction_to_be_retained
  sum_counts = 0
  for required_counts in box_counts_sorted:
    sum_counts += required_counts
    if (sum_counts > threshold):
      break
  del box_counts_sorted
  del perm
  site_selection = flex.size_t()
  for i_box in (box_counts >= required_counts).iselection():
    site_selection.extend(box_members[i_box])
  return site_selection
