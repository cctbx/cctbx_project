from __future__ import absolute_import, division, print_function
# -*- coding: utf-8 -*-
from cctbx.array_family import flex
import boost_adaptbx.boost.python as bp
from six.moves import range
from six.moves import zip
ext = bp.import_ext("cctbx_crystal_ext")
from cctbx_crystal_ext import *
from cctbx.crystal.find_best_cell import find_best_cell
from cctbx import sgtbx
from cctbx import uctbx
from cctbx import covariance, geometry
from libtbx.containers import OrderedDict
from libtbx.forward_compatibility import object
from scitbx.array_family import shared
from scitbx import stl
import scitbx.stl.set
import scitbx.stl.vector
import libtbx
from libtbx.utils import Keep, Sorry
import sys
from libtbx.utils import format_float_with_standard_uncertainty \
     as format_float_with_su
import math
from scitbx import matrix
import scitbx.cubicle_neighbors
cubicles_max_memory_allocation_set(
  number_of_bytes=scitbx.cubicle_neighbors.cubicles_max_memory_allocation_get())

pair_sym_ops = sgtbx.stl_vector_rt_mx

pair_asu_j_sym_groups = scitbx.stl.vector.set_unsigned
pair_asu_j_sym_group = scitbx.stl.set.unsigned

class symmetry(object):
  """This class represents the symmetry of a crystal and bundles information on its unit cell and space group.
  """
  def __setstate__(self, state):
    if sys.version_info.major > 2:
      from libtbx.easy_pickle import fix_py2_pickle_orig
      state = fix_py2_pickle_orig(state)
    for name,value in state.items():
      setattr(self, name, value)

  def __init__(self,
        unit_cell=None,
        space_group_symbol=None,
        space_group_info=None,
        space_group=None,
        correct_rhombohedral_setting_if_necessary=False,
        assert_is_compatible_unit_cell=True,
        raise_sorry_if_incompatible_unit_cell=False,
        force_compatible_unit_cell=True):
    """Initialises a new crystal.symmetry class object from different input data. Only one of space_group, space_group_info and space_group_symbol may be used.

    :param unit_cell:          object specifying the unit_cell properties
    :type unit_cell:           cctbc.uctbx.ext.unit_cell or tuple(lattice parameters)
    :param space_group_symbol: Hermann-Mauguin symbol of the crystallographic space group
    :type space_group_symbol:  string
    :param space_group_info:   object describing the desired space group
    :type space_group_info:    cctbx.sgtbx.space_group_info
    :param space_group:        the desired space group of the symmetry class
    :type space_group:         cctbx.sgtbx.space_group
    :param correct_rhombohedral_setting_if_necessary: If set to 'True' an automatic conversion between rhombohedral and hexagonal basis will be done
    :type correct_rhombohedral_setting_if_necessary:  boolean
    :param assert_is_compatible_unit_cell: If set to 'True' a consistency check will be performed on the relation of space group to lattice parameters
    :type assert_is_compatible_unit_cell:  boolean
    :param force_compatible_unit_cell:     If set to 'True' the crystal parameters will be averaged to comply with the restrictions of the space group
    :type force_compatible_unit_cell:      boolean

    :returns: a new crystal.symmetry class with desired properties
    :rtype: cctbx.crystal.symmetry
    """
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
      if (correct_rhombohedral_setting_if_necessary):
        sgi = self.space_group_info()
        z = sgi.group().conventional_centring_type_symbol()
        if (z == "P"):
          if (self.unit_cell().is_conventional_hexagonal_basis()):
            ls = sgi.type().lookup_symbol()
            if (ls.endswith(" :R")):
              self._space_group_info = sgi.reference_setting()
        elif (z == "R"):
          if (self.unit_cell().is_conventional_rhombohedral_basis()):
            ls = sgi.type().lookup_symbol()
            if (ls.endswith(" :H")):
              self._space_group_info = sgi.primitive_setting()
      if (assert_is_compatible_unit_cell):
        if (raise_sorry_if_incompatible_unit_cell):
          if (not self.is_compatible_unit_cell()):
            raise Sorry(("The space group '%s' is incompatible with unit cell "+
              "parameters %s.") % (str(self._space_group_info),
              str(self._unit_cell.parameters())))
        else :
          assert self.is_compatible_unit_cell(), \
            "Space group is incompatible with unit cell parameters."
      if (force_compatible_unit_cell):
        self._unit_cell = self.space_group().average_unit_cell(
          self._unit_cell)

  def _copy_constructor(self, other):
    self._unit_cell = other._unit_cell
    self._space_group_info = other._space_group_info

  def customized_copy(self, unit_cell=Keep, space_group_info=Keep,
      raise_sorry_if_incompatible_unit_cell=False):
    if (unit_cell is Keep): unit_cell = self._unit_cell
    if (space_group_info is Keep): space_group_info = self._space_group_info
    return symmetry(unit_cell=unit_cell, space_group_info=space_group_info,
      raise_sorry_if_incompatible_unit_cell=\
        raise_sorry_if_incompatible_unit_cell)

  def unit_cell(self):
    return self._unit_cell

  def space_group_number(self):
    sgi = self._space_group_info
    if sgi and sgi.group() and sgi.group().info() and \
       sgi.group().info().type() and \
       (sgi.group().info().type().number() is not None):
      return sgi.group().info().type().number()
    else:
      return None

  def space_group_info(self):
    return self._space_group_info

  def space_group(self):
    sgi = self._space_group_info
    if (sgi is None): return None
    return sgi.group()

  def as_py_code(self, indent=""):
    fmt = (
      'crystal.symmetry(\n'
      '%s  unit_cell=%s,\n'
      '%s  space_group_symbol=%s\n'
      '%s)')
    return fmt % (
      indent,
      "(%.10g, %.10g, %.10g, %.10g, %.10g, %.10g)" % self.unit_cell().parameters()
      if self.unit_cell() is not None else None,
      indent,
      '"%s"' % self.space_group_info().type().lookup_symbol()
      if self.space_group_info() is not None else None,
      indent
    )

  def show_summary(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    print(self.as_str(prefix=prefix), file=f)

  def as_str(self, prefix=""):
    return (
      prefix + "Unit cell: %s\n" % self.unit_cell() +
      prefix + "Space group: %s" % (
        self.space_group_info().symbol_and_number()
        if self.space_group_info() is not None else None)
    )

  def __str__(self):
    return self.as_str()

  def __repr__(self):
    return self.as_py_code(indent="  ")

  def is_identical_symmetry(self, other):
    ''' True if identical for self and other '''
    return self.is_similar_symmetry(other,
      relative_length_tolerance = 0,
      absolute_angle_tolerance = 0,
      absolute_length_tolerance = 0,)

  def is_similar_symmetry(self,
                          other,
                          relative_length_tolerance=0.01,
                          absolute_angle_tolerance=1.,
                          absolute_length_tolerance=-9999.,
                          ):
    if (self.unit_cell() and other.unit_cell() and
        not self.unit_cell().is_similar_to(
        other.unit_cell(),
        relative_length_tolerance,
        absolute_angle_tolerance,
        absolute_length_tolerance,
        )): return False
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
    p2n_op = p_cell.change_of_basis_op_to_niggli_cell()
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

  def join_symmetry(self, other_symmetry, force=False,
      raise_sorry_if_incompatible_unit_cell=False):
    same_result = symmetry(
        unit_cell=self.unit_cell(),
        space_group_info=self.space_group_info(),
        raise_sorry_if_incompatible_unit_cell=
          raise_sorry_if_incompatible_unit_cell)
    if (other_symmetry is None):
      return same_result
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
      space_group_info=space_group_info,
      raise_sorry_if_incompatible_unit_cell=
        raise_sorry_if_incompatible_unit_cell)

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

  def miller_set(self, indices, anomalous_flag):
    from cctbx import miller
    return miller.set(
      crystal_symmetry=self,
      indices=indices,
      anomalous_flag=anomalous_flag)

  def build_miller_set(self, anomalous_flag, d_min, d_max=None):
    from cctbx import miller
    return miller.build_set(
      crystal_symmetry=self,
      anomalous_flag=anomalous_flag,
      d_min=d_min,
      d_max=d_max)

  def as_pdb_remark_290(self):
    raise Sorry("do not use this method - not tested yet!")
    outl = ""
    if ( self.space_group() is not None and
         self.unit_cell() is not None
         ):
      uc = self.unit_cell()
      params = list(uc.parameters())
      for i,s in enumerate(self.space_group()):
        mat = s.as_double_array()
        for j in range(3):
          outl += "REMARK 290   SMTRY%d%4d%10.6f%10.6f%10.6f%15.5f\n" %(
            j+1,
            i+1,
            mat[j*3],
            mat[j*3+1],
            mat[j*3+2],
            mat[j+9]*params[j],
          )
    return outl

  def as_cif_block(self,
      cell_covariance_matrix=None,
      format="mmcif",
      numeric_format="%.3f"):
    from iotbx.cif import model
    wformat = format.lower()
    assert wformat in ("corecif", "mmcif")
    if wformat == "mmcif":
      separator = '.'
    else:
      separator = '_'
    assert numeric_format.startswith("%")
    cif_block = model.block()
    cell_prefix = '_cell%s' % separator
    if self.space_group() is not None:
      sym_loop = model.loop(data=OrderedDict((
        ('_space_group_symop'+separator+'id',
         range(1, len(self.space_group())+1)),
        ('_space_group_symop'+separator+'operation_xyz',
         [s.as_xyz() for s in self.space_group()]))))
      cif_block.add_loop(sym_loop)
      sg_prefix = '_space_group%s' % separator
      sg_type = self.space_group_info().type()
      sg = sg_type.group()
      cif_block[sg_prefix+'crystal_system'] = sg.crystal_system().lower()
      cif_block[sg_prefix+'IT_number'] = sg_type.number()
      cif_block[sg_prefix+'name_H-M_alt'] = sg_type.lookup_symbol()
      cif_block[sg_prefix+'name_Hall'] = sg_type.hall_symbol()

      sg_prefix = '_symmetry%s' % separator
      cif_block[sg_prefix+'space_group_name_H-M'] = sg_type.lookup_symbol()
      cif_block[sg_prefix+'space_group_name_Hall'] = sg_type.hall_symbol()
      cif_block[sg_prefix+'Int_Tables_number'] = sg_type.number()

    if self.unit_cell() is not None:
      uc = self.unit_cell()
      params = list(uc.parameters())
      volume = uc.volume()
      if cell_covariance_matrix is not None:
        diag = cell_covariance_matrix.matrix_packed_u_diagonal()
        for i in range(6):
          if diag[i] > 0:
            params[i] = format_float_with_su(params[i], math.sqrt(diag[i]))
        d_v_d_params = matrix.row(uc.d_volume_d_params())
        vcv = matrix.sqr(
          cell_covariance_matrix.matrix_packed_u_as_symmetric())
        var_v = (d_v_d_params * vcv).dot(d_v_d_params)
        volume = format_float_with_su(volume, math.sqrt(var_v))
        numeric_format = "%s"
      a,b,c,alpha,beta,gamma = params
      cif_block[cell_prefix+'length_a'] = numeric_format % a
      cif_block[cell_prefix+'length_b'] = numeric_format % b
      cif_block[cell_prefix+'length_c'] = numeric_format % c
      cif_block[cell_prefix+'angle_alpha'] = numeric_format % alpha
      cif_block[cell_prefix+'angle_beta'] = numeric_format % beta
      cif_block[cell_prefix+'angle_gamma'] = numeric_format % gamma
      cif_block[cell_prefix+'volume'] = numeric_format % volume
    return cif_block

  def expand_to_p1(self, sites_cart):
    from scitbx import matrix
    result = flex.vec3_double()
    for site_cart in sites_cart:
      for smx in self.space_group().smx():
        m3 = smx.r().as_double()
        m3 = matrix.sqr(m3)
        t = smx.t().as_double()
        t = matrix.col((t[0],t[1],t[2]))
        site_frac=flex.vec3_double([self.unit_cell().fractionalize(site_cart),])
        result.append(self.unit_cell().orthogonalize(m3.elems*site_frac+t)[0])
    return result

  def is_nonsense(self):
    uc = self.unit_cell()
    if uc is None:
      return True
    if(self.space_group_info() is None or self.space_group() is None):
      return True
    ucp = uc.parameters()
    result = ((abs(1.-ucp[0])<1.e-3 and
             abs(1.-ucp[1])<1.e-3 and
             abs(1.-ucp[2])<1.e-3) or
          (abs(0.-ucp[0])<1.e-3 and
             abs(0.-ucp[1])<1.e-3 and
             abs(0.-ucp[2])<1.e-3) or
          (abs(1.-ucp[3])<1.e-3 and
             abs(1.-ucp[4])<1.e-3 and
             abs(1.-ucp[5])<1.e-3) or
          (abs(0.-ucp[3])<1.e-3 and
             abs(0.-ucp[4])<1.e-3 and
             abs(0.-ucp[5])<1.e-3))
    return result

  def is_empty(self):
    return self.unit_cell() is None and self.space_group_info() is None

  def is_incomplete(self):
    return self.unit_cell() is None or self.space_group() is None

def select_crystal_symmetry(
      from_command_line     = None,
      from_parameter_file   = None,
      from_coordinate_files = [None],
      from_reflection_files = [None],
      enforce_similarity = False,
      absolute_angle_tolerance = 1.,
      absolute_length_tolerance = -9999.):
  """Select/construct a crystal symmetry from a list of various options"""
  tmp = [from_command_line, from_parameter_file]+from_coordinate_files \
        +from_reflection_files
  if tmp.count(None)==len(tmp):
    raise AssertionError("No unit cell and symmetry information supplied")
  if len(tmp) > 0 and enforce_similarity:
    cs0 = None
    i = 0
    while cs0 is None and i<len(tmp):
      cs0 = tmp[i]
      if cs0 is not None:
        if cs0.is_nonsense() or cs0.is_empty():
          cs0 = None
      i += 1
    for cs in tmp[i:]:
      if cs and not cs.is_nonsense() and not cs.is_empty():
        is_similar_cs = cs0.is_similar_symmetry(cs,
           absolute_angle_tolerance=absolute_angle_tolerance,
           absolute_length_tolerance=absolute_length_tolerance)
        if(not is_similar_cs):
          msg = "Crystal symmetry mismatch between different files.\n"
          msg += "%s %s\n" % (cs0.unit_cell(), cs0.space_group_info())
          msg += "%s %s\n" % (cs.unit_cell(), cs.space_group_info())
          raise Sorry("%s"%(msg))
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
      if enforce_similarity:  # usual, require compatibility here
        result = result.join_symmetry(
          other_symmetry=crystal_symmetry, force=False)
      else:  # skip incompatible symmetries (can happen e.g. in map_box)
        try:
          result = result.join_symmetry(
            other_symmetry=crystal_symmetry, force=False)
        except Exception as e:
          pass
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
  if result.is_nonsense():
    return None
  if result.is_empty():
    return None
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

  def xray_structure(self, scatterers=None):
    import cctbx.xray.structure as _
    return _(special_position_settings=self, scatterers=scatterers)

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
      error_message="Excessive special position correction:"):
  assert (site_frac is None) != (site_cart is None)
  unit_cell = crystal_symmetry.unit_cell()
  if (site_frac is None):
    site_frac = unit_cell.fractionalize(site_cart)
  site_special_frac = special_op * site_frac
  distance_moved = unit_cell.distance(site_special_frac, site_frac)
  if (distance_moved > tolerance):
    error_message += "\n  unit_cell: %s" % str(unit_cell)
    error_message += "\n  space_group_info: %s" % str(
      crystal_symmetry.space_group_info())
    error_message += "\n  special_op: %s" % str(special_op)
    if (site_label is not None):
      error_message += "\n  site_label: %s" % site_label
    error_message += "\n  site_frac: %s" % str(site_frac)
    error_message += "\n  site_special_frac: %s" % str(site_special_frac)
    error_message += "\n  distance_moved: %g" % distance_moved
    raise AssertionError(error_message)
  if (site_cart is None):
    return site_special_frac
  return unit_cell.orthogonalize(site_special_frac)

@bp.inject_into(pair_asu_table)
class _():

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
        print("i_seq:", i_seq, file=f)
        for j_seq,j_sym_group in j_seq_dict.items():
          print("  j_seq:", j_seq, file=f)
          for j_syms in j_sym_group:
            print("    j_syms:", list(j_syms), file=f)
    else:
      assert len(site_labels) == self.table().size()
      for i_seq, j_seq_dict in enumerate(self.table()):
        print("%s(%d)" % (site_labels[i_seq], i_seq), file=f)
        for j_seq,j_sym_group in j_seq_dict.items():
          print("  %s(%d)" % (site_labels[j_seq], j_seq), file=f)
          for j_syms in j_sym_group:
            print("    j_syms:", list(j_syms), file=f)

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

  def show_dihedral_angles(self,
        site_labels=None,
        sites_frac=None,
        sites_cart=None,
        keep_pair_asu_table=False,
        max_d=1.7,
        max_angle=170,
        out=None):
    return show_dihedral_angles(
      pair_asu_table=self,
      site_labels=site_labels,
      sites_frac=sites_frac,
      sites_cart=sites_cart,
      max_d=max_d,
      max_angle=max_angle,
      out=out)

class calculate_distances(object):

  def __init__(self,
               pair_asu_table,
               sites_frac,
               skip_j_seq_less_than_i_seq=True,
               covariance_matrix=None,
               cell_covariance_matrix=None,
               parameter_map=None):
    libtbx.adopt_init_args(self, locals())
    self.distances = flex.double()
    if self.covariance_matrix is not None:
      self.variances = flex.double()
    else:
      self.variances = None
    self.pair_counts = flex.size_t()

  def __iter__(self):
    return next(self)

  def __next__(self):

    class distance(object):
      def __init__(self,
                   distance,
                   i_seq,
                   j_seq,
                   pair_count,
                   rt_mx_ji=None,
                   i_j_sym=None,
                   variance=None):
        libtbx.adopt_init_args(self, locals())

    asu_mappings = self.pair_asu_table.asu_mappings()
    unit_cell = asu_mappings.unit_cell()
    if self.covariance_matrix is not None:
      assert self.parameter_map is not None
      cov_cart = covariance.orthogonalize_covariance_matrix(
        self.covariance_matrix, unit_cell, self.parameter_map)
    for i_seq,asu_dict in enumerate(self.pair_asu_table.table()):
      rt_mx_i_inv = asu_mappings.get_rt_mx(i_seq, 0).inverse()
      site_frac_i = self.sites_frac[i_seq]
      pair_count = 0
      dists = flex.double()
      j_seq_i_group = []
      for j_seq,j_sym_groups in asu_dict.items():
        if self.skip_j_seq_less_than_i_seq and j_seq < i_seq: continue
        site_frac_j = self.sites_frac[j_seq]
        for i_group,j_sym_group in enumerate(j_sym_groups):
          pair_count += j_sym_group.size()
          j_sym = j_sym_group[0]
          rt_mx_ji = rt_mx_i_inv.multiply(asu_mappings.get_rt_mx(j_seq, j_sym))
          dist = unit_cell.distance(site_frac_i, rt_mx_ji * site_frac_j)
          dists.append(dist)
          j_seq_i_group.append((j_seq,i_group))
      permutation = flex.sort_permutation(data=dists)
      for j_seq,i_group in flex.select(j_seq_i_group, permutation):
        site_frac_j = self.sites_frac[j_seq]
        j_sym_groups = asu_dict[j_seq]
        j_sym_group = j_sym_groups[i_group]
        for i_j_sym,j_sym in enumerate(j_sym_group):
          rt_mx_ji = rt_mx_i_inv.multiply(
            asu_mappings.get_rt_mx(j_seq, j_sym))
          site_frac_ji = rt_mx_ji * site_frac_j
          d = geometry.distance((unit_cell.orthogonalize(site_frac_i),
                                 unit_cell.orthogonalize(site_frac_ji)))
          dist = d.distance_model
          self.distances.append(dist)
          if self.covariance_matrix is not None:
            cov = covariance.extract_covariance_matrix_for_sites(
              flex.size_t((i_seq,j_seq)), cov_cart, self.parameter_map)
            if self.cell_covariance_matrix is not None:
              var = d.variance(cov, unit_cell, rt_mx_ji)
              # a bit of a hack - if the variance is approx zero, we assume
              # that the distance is constrained (e.g. as part of a rigid body)
              # unless it is a distance between two symmetry-related atoms, in
              # which case we include the cell errors into the variance calculation
              if var > 2e-16 or (i_seq == j_seq and not rt_mx_ji.is_unit_mx()):
                var = d.variance(
                  cov, self.cell_covariance_matrix, unit_cell, rt_mx_ji)
            else:
              var = d.variance(cov, unit_cell, rt_mx_ji)
            self.variances.append(var)
          else:
            var = None
          yield distance(
            dist, i_seq, j_seq, pair_count, rt_mx_ji, i_j_sym, variance=var)

      self.pair_counts.append(pair_count)

class show_distances(libtbx.slots_getstate_setstate):

  __slots__ = ["pair_asu_table", "distances_info", "have_sym"]

  def __init__(self,
        pair_asu_table,
        site_labels=None,
        sites_frac=None,
        sites_cart=None,
        show_cartesian=False,
        keep_pair_asu_table=False,
        skip_j_seq_less_than_i_seq=False,
        out=None):
    assert [sites_frac, sites_cart].count(None) == 1
    if (out is None): out = sys.stdout
    if (keep_pair_asu_table):
      self.pair_asu_table = pair_asu_table
    else:
      self.pair_asu_table = None
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
    self.distances_info = calculate_distances(
      pair_asu_table, sites_frac,
      skip_j_seq_less_than_i_seq=skip_j_seq_less_than_i_seq)
    self.have_sym = False
    i_seqs_done = set()
    for di in self.distances_info:
      i_seq, j_seq = di.i_seq, di.j_seq
      rt_mx_ji = di.rt_mx_ji
      from scitbx.array_family import flex
      first_time_i_seq = (i_seq not in i_seqs_done)
      if (first_time_i_seq):
        i_seqs_done.add(i_seq)
        if (site_labels is None):
          s = label_fmt % (i_seq+1)
        else:
          s = label_fmt % site_labels[i_seq]
        s += " pair count: %3d" % di.pair_count
        site_frac_i = sites_frac[i_seq]
        if (show_cartesian):
          formatted_site = [" %7.2f" % x
            for x in unit_cell.orthogonalize(site_frac_i)]
        else:
          formatted_site = [" %7.4f" % x for x in site_frac_i]
        print(("%%-%ds" % (label_len+23)) % s, \
          "<<"+",".join(formatted_site)+">>", file=out)
      if (site_labels is None):
        print(" ", label_fmt % (j_seq+1) + ":", end=' ', file=out)
      else:
        print(" ", label_fmt % (site_labels[j_seq] + ":"), end=' ', file=out)
      print("%8.4f" % di.distance, end=' ', file=out)
      if di.i_j_sym != 0:
        s = "sym. equiv."
      else:
        s = "           "
      site_frac_ji = rt_mx_ji * sites_frac[j_seq]
      if (show_cartesian):
        formatted_site = [" %7.2f" % x
          for x in unit_cell.orthogonalize(site_frac_ji)]
      else:
        formatted_site = [" %7.4f" % x for x in site_frac_ji]
      s += " (" + ",".join(formatted_site) +")"
      if (not rt_mx_ji.is_unit_mx()):
        s += " sym=" + str(rt_mx_ji)
        self.have_sym = True
      print(s, file=out)
      if first_time_i_seq and di.pair_count == 0:
        print("  no neighbors", file=out)

class calculate_angles(object):

  def __init__(self,
               pair_asu_table,
               sites_frac,
               skip_j_seq_less_than_i_seq=True,
               covariance_matrix=None,
               cell_covariance_matrix=None,
               parameter_map=None,
               conformer_indices=None):
    libtbx.adopt_init_args(self, locals())
    self.distances = flex.double()
    if self.covariance_matrix is not None:
      self.variances = flex.double()
    else:
      self.variances = None
    self.angles = flex.double()
    self.pair_counts = flex.size_t()

  def __iter__(self):
    return next(self)

  def __next__(self):

    class angle(object):
      def __init__(self,
                   angle,
                   i_seqs,
                   rt_mx_ji=None,
                   rt_mx_ki=None,
                   variance=None):
        libtbx.adopt_init_args(self, locals())

    asu_mappings = self.pair_asu_table.asu_mappings()
    unit_cell = asu_mappings.unit_cell()
    if self.covariance_matrix is not None:
      assert self.parameter_map is not None
      cov_cart = covariance.orthogonalize_covariance_matrix(
        self.covariance_matrix, unit_cell, self.parameter_map)

    ## angle is formed by j_seq-i_seq-k_seq
    for i_seq,asu_dict in enumerate(self.pair_asu_table.table()):
      rt_mx_i_inv = asu_mappings.get_rt_mx(i_seq, 0).inverse()
      site_frac_i = self.sites_frac[i_seq]
      angles = flex.double()
      for j_seq,j_sym_groups in asu_dict.items():
        site_frac_j = self.sites_frac[j_seq]
        for j_sym_group in j_sym_groups:
          for i_j_sym,j_sym in enumerate(j_sym_group):
            rt_mx_ji = rt_mx_i_inv.multiply(
              asu_mappings.get_rt_mx(j_seq, j_sym))
            site_frac_ji = rt_mx_ji * site_frac_j
            for k_seq, k_sym_groups in asu_dict.items():
              if self.skip_j_seq_less_than_i_seq and j_seq < k_seq: continue
              if k_seq == j_seq and j_sym_group.size() <= 1: continue
              if k_seq > j_seq: continue
              site_frac_k = self.sites_frac[k_seq]
              for k_sym_group in k_sym_groups:
                for i_k_sym,k_sym in enumerate(k_sym_group):
                  if j_seq == k_seq and i_j_sym <= i_k_sym: continue
                  if i_seq == k_seq and i_k_sym == 0: continue
                  if (self.conformer_indices is not None and
                      self.conformer_indices[j_seq] !=
                      self.conformer_indices[k_seq]): continue
                  rt_mx_ki = rt_mx_i_inv.multiply(
                    asu_mappings.get_rt_mx(k_seq, k_sym))
                  site_frac_ki = rt_mx_ki * site_frac_k
                  a = geometry.angle((unit_cell.orthogonalize(site_frac_ji),
                                      unit_cell.orthogonalize(site_frac_i),
                                      unit_cell.orthogonalize(site_frac_ki)))
                  angle_ = a.angle_model
                  self.angles.append(angle_)
                  if self.covariance_matrix is not None:
                    cov = covariance.extract_covariance_matrix_for_sites(
                      flex.size_t((j_seq, i_seq, k_seq)),
                      cov_cart, self.parameter_map)
                    if self.cell_covariance_matrix is not None:
                      var = a.variance(
                        cov, unit_cell, (rt_mx_ji, sgtbx.rt_mx(), rt_mx_ki))
                      # a bit of a hack - see equivalent comment in
                      # calculate_distances
                      if (var > 2e-15 or
                          (i_seq == j_seq and not rt_mx_ji.is_unit_mx()) or
                          (i_seq == k_seq and not rt_mx_ki.is_unit_mx())):
                        var = a.variance(
                          cov, self.cell_covariance_matrix, unit_cell,
                          (rt_mx_ji, sgtbx.rt_mx(), rt_mx_ki))
                    else:
                      var = a.variance(
                        cov, unit_cell, (rt_mx_ji, sgtbx.rt_mx(), rt_mx_ki))
                    self.variances.append(var)
                  else:
                    var = None
                  yield angle(angle_, (j_seq, i_seq, k_seq),
                              rt_mx_ji, rt_mx_ki, variance=var)

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
    rt_mxs = []
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
    angles = calculate_angles(pair_asu_table, sites_frac)
    for a in angles:
      j_seq, i_seq, k_seq = a.i_seqs
      rt_mx_ji = a.rt_mx_ji
      rt_mx_ki = a.rt_mx_ki
      if (site_labels is None):
        s = label_fmt % (j_seq+1) + ":"
      else:
        i_label = site_labels[i_seq]
        j_label = site_labels[j_seq]
        k_label = site_labels[k_seq]
        if not rt_mx_ji.is_unit_mx():
          if rt_mx_ji in rt_mxs:
            j = rt_mxs.index(rt_mx_ji) + 1
          else:
            rt_mxs.append(rt_mx_ji)
            j = len(rt_mxs)
          j_label += "*%s" %j
        if not rt_mx_ki.is_unit_mx():
          if rt_mx_ki in rt_mxs:
            k = rt_mxs.index(rt_mx_ki) + 1
          else:
            rt_mxs.append(rt_mx_ki)
            k = len(rt_mxs)
          k_label += "*%s" %k
        s = label_fmt % (j_label, i_label, k_label)
      s += " %6.2f" % a.angle
      print(s, file=out)

    self.angles = angles.angles
    self.distance = angles.distances
    for i, rt_mx in enumerate(rt_mxs):
      print("*%s" %(i+1), end=' ', file=out)
      print(rt_mx, file=out)


class dihedral_angle_def(object):
  def __init__(self, seqs, rt_mxs):
    libtbx.adopt_init_args(self, locals())

class calculate_dihedrals(object):
  def __init__(self,
               pair_asu_table,
               sites_frac,
               dihedral_defs=None, #angle definition
               skip_j_seq_less_than_i_seq=True,
               covariance_matrix=None,
               cell_covariance_matrix=None,
               parameter_map=None,
               conformer_indices=None,
               max_d = 1.9,
               max_angle = 170):
    libtbx.adopt_init_args(self, locals())
    if self.covariance_matrix is not None:
      self.variances = flex.double()
    else:
      self.variances = None
    self.dihedrals = flex.double()

  def __iter__(self):
    return self.next()

  def next(self):

    class dihedral(object):
      def __init__(self,
                   angle,
                   i_seqs,
                   rt_mxs,
                   variance=None):
        libtbx.adopt_init_args(self, locals())
      def __eq__(self, other):
        for i in range(0,4):
          if self.i_seqs[i] != other.i_seqs[i] or self.rt_mxs[i] != other.rt_mxs[i]:
            return False
        return True

    asu_mappings = self.pair_asu_table.asu_mappings()
    unit_cell = asu_mappings.unit_cell()
    if self.covariance_matrix is not None:
      assert self.parameter_map is not None
      cov_cart = covariance.orthogonalize_covariance_matrix(
        self.covariance_matrix, unit_cell, self.parameter_map)
    if self.dihedral_defs is not None:
      for ad in self.dihedral_defs:
        sites = []
        for i in range(0,4):
          site_frac = ad.rt_mxs[i] * self.sites_frac[ad.seqs[i]]
          sites.append(unit_cell.orthogonalize(site_frac))
        a = geometry.dihedral(sites)
        angle_ = a.dihedral_model
        self.dihedrals.append(angle_)
        if self.covariance_matrix is not None:
          cov = covariance.extract_covariance_matrix_for_sites(
            flex.size_t(ad.seqs), cov_cart, self.parameter_map)
          if self.cell_covariance_matrix is not None:
            var = a.variance(cov, unit_cell, ad.rt_mxs)
          else:
            var = a.variance(cov, unit_cell, ad.rt_mxs)
          self.variances.append(var)
        else:
          var = None
        yield dihedral(angle_, ad.seqs, ad.rt_mxs, variance=var)
      return
    table = self.pair_asu_table.table()
    for i_seq,i_asu_dict in enumerate(table):
      rt_mx_i_inv = asu_mappings.get_rt_mx(i_seq, 0).inverse()
      i_site_frac = self.sites_frac[i_seq]
      for j_seq,j_sym_groups in i_asu_dict.items():
        if j_seq < i_seq: continue
        rt_mx_j0_inv = asu_mappings.get_rt_mx(j_seq, 0).inverse()
        for j_sym_group in j_sym_groups:
          for j_sym_idx,j_sym in enumerate(j_sym_group):
            rt_mx_j = rt_mx_i_inv.multiply(asu_mappings.get_rt_mx(j_seq, j_sym))
            j_site_frac = rt_mx_j * self.sites_frac[j_seq]
            for k_seq, k_sym_groups in table[j_seq].items():
              if k_seq < i_seq: continue
              if (self.conformer_indices is not None and
                  self.conformer_indices[j_seq] !=
                  self.conformer_indices[k_seq]):
                continue
              rt_mx_k0_inv = asu_mappings.get_rt_mx(k_seq, 0).inverse()
              rt_mx_kj = rt_mx_j.multiply(rt_mx_j0_inv)
              for k_sym_group in k_sym_groups:
                for k_sym_idx,k_sym in enumerate(k_sym_group):
                  rt_mx_k = rt_mx_kj.multiply(asu_mappings.get_rt_mx(k_seq, k_sym))
                  k_site_frac = rt_mx_k *  self.sites_frac[k_seq]
                  if k_site_frac == i_site_frac:
                    continue
                  if unit_cell.distance(j_site_frac, k_site_frac) > self.max_d:
                    continue
                  rt_mx_lk = rt_mx_k.multiply(rt_mx_k0_inv)
                  for l_seq, l_sym_groups in table[k_seq].items():
                    if l_seq < i_seq: continue
                    if (self.conformer_indices is not None and
                        self.conformer_indices[k_seq] !=
                        self.conformer_indices[l_seq]):
                      continue
                    for l_sym_group in l_sym_groups:
                      for l_sym_idx, l_sym in enumerate(l_sym_group):
                        rt_mx_l = rt_mx_lk.multiply(asu_mappings.get_rt_mx(l_seq, l_sym))
                        l_site_frac = rt_mx_l *  self.sites_frac[l_seq]
                        if l_site_frac in (i_site_frac, j_site_frac):
                          continue
                        if unit_cell.angle(i_site_frac, j_site_frac, k_site_frac) > self.max_angle or\
                           unit_cell.angle(j_site_frac, k_site_frac, l_site_frac) > self.max_angle:
                          continue
                        rt_mxs = [sgtbx.rt_mx(), rt_mx_j, rt_mx_k, rt_mx_l]
                        sites = [i_site_frac, j_site_frac, k_site_frac, l_site_frac]
                        seqs = [i_seq, j_seq, k_seq, l_seq]
                        rtmx_count = [0,0,0,0]
                        # find the most common matrix to eliminate
                        for idx, rt_mx in enumerate(rt_mxs):
                          for idx1 in range(idx+1, 4):
                            if rt_mx == rt_mxs[idx1]:
                              rtmx_count[idx] += 1
                        max_idx = rtmx_count.index(max(rtmx_count))
                        rt_mx_inv = rt_mxs[max_idx].inverse()
                        for i in range(0, 4):
                          sites[i] = rt_mxs[i].inverse() * sites[i]
                          rt_mxs[i] = rt_mx_inv.multiply(rt_mxs[i])
                          sites[i] = unit_cell.orthogonalize(rt_mxs[i] * sites[i])
                        a = geometry.dihedral(sites)
                        angle_ = a.dihedral_model
                        self.dihedrals.append(angle_)
                        if self.covariance_matrix is not None:
                          cov = covariance.extract_covariance_matrix_for_sites(
                            flex.size_t(seqs), cov_cart, self.parameter_map)
                          if self.cell_covariance_matrix is not None:
                            var = a.variance(cov, unit_cell, rt_mxs)
                          else:
                            var = a.variance(cov, unit_cell, rt_mxs)
                          self.variances.append(var)
                        else:
                          var = None
                        yield dihedral(angle_, seqs, rt_mxs, variance=var)

class show_dihedral_angles(object):

  def __init__(self,
        pair_asu_table,
        site_labels=None,
        sites_frac=None,
        sites_cart=None,
        show_cartesian=False,
        max_d=1.7,
        max_angle=170,
        out=None):

    assert [sites_frac, sites_cart].count(None) == 1
    if (out is None): out = sys.stdout
    rt_mxs = []
    if (site_labels is None):
      label_len = len("%d" % (sites_frac.size()+1))
      label_fmt = "site_%%0%dd" % label_len
      label_len += 5
    else:
      label_len = 1
      for label in site_labels:
        label_len = max(label_len, len(label))
      label_fmt = "%%-%ds" % (label_len+4)
    angles = calculate_dihedrals(pair_asu_table, sites_frac,
      max_d=max_d, max_angle=max_angle)
    for d in angles:
      if (site_labels is None):
        s = label_fmt % (d.i_seqs[0]+1) + ":"
      else:
        s = ""
        for idx, i_seq in enumerate(d.i_seqs):
          label = label_fmt % site_labels[i_seq]
          rt_mx = d.rt_mxs[idx]
          if not rt_mx.is_unit_mx():
            if rt_mx in rt_mxs:
              j = rt_mxs.index(rt_mx) + 1
            else:
              rt_mxs.append(rt_mx)
              j = len(rt_mxs)
            label += "*%s" %j
          s += label
      s += " %6.2f" % d.angle
      print(s, file=out)

    self.dihedrals = angles.dihedrals
    for i, rt_mx in enumerate(rt_mxs):
      print("*%s %s" %(i+1, rt_mx), file=out)

class sym_pair(libtbx.slots_getstate_setstate):

  __slots__ = ["i_seq", "j_seq", "rt_mx_ji"]

  def __init__(self, i_seq, j_seq, rt_mx_ji):
    self.i_seq = i_seq
    self.j_seq = j_seq
    self.rt_mx_ji = rt_mx_ji

  def i_seqs(self):
    return (self.i_seq, self.j_seq)

@bp.inject_into(pair_sym_table)
class _():

  def iterator(self):
    for i_seq,pair_sym_dict in enumerate(self):
      for j_seq,sym_ops in pair_sym_dict.items():
        for rt_mx_ji in sym_ops:
          yield sym_pair(i_seq=i_seq, j_seq=j_seq, rt_mx_ji=rt_mx_ji)

  def tidy(self, site_symmetry_table):
    result = pair_sym_table(size=self.size())
    for i_seq,pair_sym_dict in enumerate(self):
      for j_seq,sym_ops in pair_sym_dict.items():
        sepi_objs = []
        for rt_mx_ji in sym_ops:
          if (i_seq <= j_seq):
            i, j = i_seq, j_seq
          else:
            i, j, rt_mx_ji = j_seq, i_seq, rt_mx_ji.inverse()
          for sepi_obj in sepi_objs:
            if (sepi_obj.is_equivalent(rt_mx_ji=rt_mx_ji)):
              break
          else:
            sepi_obj = site_symmetry_table \
              .symmetry_equivalent_pair_interactions(
                i_seq=i, j_seq=j, rt_mx_ji=rt_mx_ji)
            sepi_objs.append(sepi_obj)
        ri = result[i]
        ri[j] = pair_sym_ops()
        rij = ri[j]
        for rt_mx_ji in sorted([sepi_obj.get()[0] for sepi_obj in sepi_objs]):
          rij.append(rt_mx_ji)
    return result

  def full_connectivity(self, site_symmetry_table=None):
    result = pair_sym_table(size=self.size())
    for i_seq,pair_sym_dict in enumerate(self):
      ri = result[i_seq]
      for j_seq,sym_ops in pair_sym_dict.items():
        rij = ri.get(j_seq)
        if (rij is None):
          ri[j_seq] = pair_sym_ops()
          rij = ri[j_seq]
        rj = result[j_seq]
        rji = rj.get(i_seq)
        if (rji is None):
          rj[i_seq] = pair_sym_ops()
          rji = rj[i_seq]
        for rt_mx_ji in sym_ops:
          if (site_symmetry_table is None):
            rij.append(rt_mx_ji)
            if (i_seq != j_seq):
              rji.append(rt_mx_ji.inverse())
          else:
            sepi = site_symmetry_table.symmetry_equivalent_pair_interactions
            for s in sepi(
                       i_seq=i_seq,
                       j_seq=j_seq,
                       rt_mx_ji=rt_mx_ji).get():
              rij.append(s)
            if (i_seq != j_seq):
              for s in sepi(
                         i_seq=j_seq,
                         j_seq=i_seq,
                         rt_mx_ji=rt_mx_ji.inverse()).get():
                rji.append(s)
    return result

  def show(self, f=None,
        site_labels=None,
        site_symmetry_table=None,
        sites_frac=None,
        unit_cell=None):
    if (site_labels is not None):
      assert len(site_labels) == self.size()
    if (sites_frac is not None):
      assert len(sites_frac) == self.size()
      assert unit_cell is not None
    if (site_symmetry_table is not None):
      assert site_symmetry_table.indices().size() == self.size()
    if (f is None): f = sys.stdout
    def show_i():
      if (site_labels is None):
        print("i_seq:", i_seq, file=f)
      else:
        print("%s(%d)" % (site_labels[i_seq], i_seq), file=f)
    def show_j():
      if (site_labels is None):
        print("  j_seq:", j_seq, file=f)
      else:
        print("  %s(%d)" % (site_labels[j_seq], j_seq), file=f)
    for i_seq,pair_sym_dict in enumerate(self):
      show_i()
      for j_seq,sym_ops in pair_sym_dict.items():
        show_j()
        if (site_symmetry_table is None):
          if (sites_frac is None):
            for rt_mx_ji in sym_ops:
              print("   ", rt_mx_ji, file=f)
          elif (len(sym_ops) > 0):
            max_len = max([len(str(_)) for _ in sym_ops])
            fmt = "    %%-%ds  %%8.4f" % max_len
            for rt_mx_ji in sym_ops:
              d = unit_cell.distance(
                sites_frac[i_seq], rt_mx_ji * sites_frac[j_seq])
              print(fmt % (str(rt_mx_ji), d), file=f)
        else:
          max_len = 0
          sepis = []
          for rt_mx_ji in sym_ops:
            sepi = site_symmetry_table.symmetry_equivalent_pair_interactions(
              i_seq=i_seq, j_seq=j_seq, rt_mx_ji=rt_mx_ji).get()
            sepis.append(sepi)
            max_len = max(max_len, max([len(str(_)) for _ in sepi]))
          if (max_len != 0):
            fmt = "    %%-%ds%%s%%s" % max_len
            for sepi in sepis:
              d = ""
              e = ""
              for s in sepi:
                if (sites_frac is not None):
                  d = "  %8.4f" % unit_cell.distance(
                    sites_frac[i_seq], s * sites_frac[j_seq])
                print((fmt % (str(s), d, e)).rstrip(), file=f)
                e = "  sym. equiv."

  def show_distances(self,
        unit_cell,
        site_symmetry_table,
        site_labels=None,
        sites_frac=None,
        sites_cart=None,
        show_cartesian=False,
        skip_j_seq_less_than_i_seq=False,
        skip_sym_equiv=False,
        out=None):
    assert [sites_frac, sites_cart].count(None) == 1
    if (out is None): out = sys.stdout
    if (sites_frac is None):
      sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
    if (site_labels is None):
      label_len = len("%d" % (sites_frac.size()+1))
      label_fmt = "site_%%0%dd" % label_len
      label_len += 5
    else:
      label_len = max(1, max([len(_) for _ in site_labels]))
      label_fmt = "%%-%ds" % (label_len+1)
    if (skip_j_seq_less_than_i_seq):
      back_interactions = None
    else:
      back_interactions = [set() for _ in range(self.size())]
      for i_seq,pair_sym_dict in enumerate(self):
        for j_seq,sym_ops in pair_sym_dict.items():
          if (j_seq != i_seq):
            assert i_seq < j_seq
            back_interactions[j_seq].add(i_seq)
    pair_counts = flex.size_t()
    for i_seq,pair_sym_dict in enumerate(self):
      site_frac_i = sites_frac[i_seq]
      distance_info = []
      def distance_info_append(j_seq, rt_mx_ji):
        site_frac_ji = rt_mx_ji * sites_frac[j_seq]
        sepi = site_symmetry_table.symmetry_equivalent_pair_interactions(
          i_seq=i_seq, j_seq=j_seq, rt_mx_ji=rt_mx_ji).get()
        distance_info.append([
          unit_cell.distance(site_frac_i, site_frac_ji),
          j_seq,
          rt_mx_ji,
          sepi])
      for j_seq,sym_ops in pair_sym_dict.items():
        for rt_mx_ji in sym_ops:
          distance_info_append(j_seq, rt_mx_ji)
      if (back_interactions is not None):
        for j_seq in back_interactions[i_seq]:
          for rt_mx_ji_inv in self[j_seq][i_seq]:
            distance_info_append(j_seq, rt_mx_ji_inv.inverse())
      distance_info.sort()
      if (skip_sym_equiv):
        pair_count = len(distance_info)
      else:
        pair_count = sum([len(sepi) for _,_,_,sepi in distance_info])
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
      print(("%%-%ds" % (label_len+23)) % s, \
        "<<"+",".join(formatted_site)+">>", file=out)
      for distance,j_seq,_,sepi in distance_info:
        sym_equiv = "           "
        for rt_mx_ji_equiv in sepi:
          if (site_labels is None):
            print(" ", label_fmt % (j_seq+1) + ":", end=' ', file=out)
          else:
            print(" ", label_fmt % (site_labels[j_seq] + ":"), end=' ', file=out)
          print("%8.4f" % distance, end=' ', file=out)
          s = sym_equiv
          sym_equiv = "sym. equiv."
          site_frac_ji_equiv = rt_mx_ji_equiv * sites_frac[j_seq]
          if (show_cartesian):
            formatted_site = [" %7.2f" % x
              for x in unit_cell.orthogonalize(site_frac_ji_equiv)]
          else:
            formatted_site = [" %7.4f" % x for x in site_frac_ji_equiv]
          s += " (" + ",".join(formatted_site) +")"
          if (not rt_mx_ji_equiv.is_unit_mx()):
            s += " sym=" + str(rt_mx_ji_equiv)
          print(s, file=out)
          if (skip_sym_equiv):
            break
      if (pair_count == 0):
        print("  no neighbors", file=out)
      pair_counts.append(pair_count)
    return pair_counts

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

  def symmetry_edge_list(self):
    result = []
    for i_seq,pair_sym_dict in enumerate(self):
      for j_seq,sym_ops in pair_sym_dict.items():
        for sym_op in sym_ops:
          if not (sym_op.is_unit_mx()):
            result.append((i_seq, j_seq, sym_op))
    return result

  def both_edge_list(self):
    simple = []
    sym = []
    for i_seq,pair_sym_dict in enumerate(self):
      for j_seq,sym_ops in pair_sym_dict.items():
        for sym_op in sym_ops:
          if (sym_op.is_unit_mx()):
            simple.append((i_seq, j_seq))
          else:
            sym.append((i_seq, j_seq, sym_op))
    return simple, sym

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

  def discard_symmetry(self):
    result = pair_sym_table()
    sym_ops = sgtbx.space_group().all_ops()
    for i_seq,self_pair_sym_dict in enumerate(self):
      d = pair_sym_dict()
      for j_seq in self_pair_sym_dict.keys():
        d[j_seq] = sym_ops
      result.append(d)
    return result

  def add_pair_sym_table_in_place(self, other):
    self_size = self.size()
    assert other.size() <= self_size
    for self_pair_sym_dict,other_pair_sym_dict in zip(self, other):
      for j_seq,other_sym_ops in other_pair_sym_dict.items():
        assert j_seq < self_size
        self_pair_sym_dict.setdefault(j_seq)
        self_pair_sym_dict[j_seq].extend(other_sym_ops)

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
    assignments = flex.size_t(range(n_sites))
    n_clusters = n_sites
    index_groups = [[i_seq] for i_seq in range(n_sites)]
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
    assignments = flex.size_t(range(n_sites))
    n_clusters = n_sites
    index_groups = [[i_seq] for i_seq in range(n_sites)]
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
  box_members = [flex.size_t() for i in range(box_counts.size())]
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

def unit_crystal_symmetry():
  # returns crystal symmetry with unit cell=(1,1,1,90,90,90), sg=1
  sg=sgtbx.space_group_info(symbol='P1')
  uc=uctbx.unit_cell((1,1,1,90,90,90))
  from cctbx import crystal
  return crystal.symmetry(unit_cell=uc,space_group_info=sg)

bp.inject(ext.neighbors_simple_pair_generator, bp.py3_make_iterator)
bp.inject_into(ext.neighbors_fast_pair_generator, bp.py3_make_iterator)
