from __future__ import absolute_import, division, print_function
import cctbx.array_family.flex # import dependency
import scitbx.array_family.shared # import dependency

import boost_adaptbx.boost.python as bp
from six.moves import zip
ext = bp.import_ext("cctbx_uctbx_ext")
from cctbx_uctbx_ext import *

from scitbx import matrix
import sys

class unit_cell(ext.unit_cell):
  """
  Class for the handling of unit cell information.

  All angles are in degrees.

  The PDB convention for orthogonalization and fractionalization
  of coordinates is used:

  | Crystallographic Basis: D = {a,b,c}
  | Cartesian Basis:        C = {i,j,k}
  | i parallel to a
  | j is in (a,b) plane
  | k = i x j

  :param parameters: A list, tuple or string of unit cell parameters.
  :param metrical_matrix: Metrical matrix. See also
                          :py:meth:`metrical_matrix`.
  :param orthogonalization_matrix: Orthogonalization matrix. See also
                          :py:meth:`orthogonalization_matrix`.

  :returns: None
  :rtype: None
  """

  def __init__(self,
        parameters=None,
        metrical_matrix=None,
        orthogonalization_matrix=None):
    assert [parameters, metrical_matrix, orthogonalization_matrix].count(None) >= 2
    if (parameters is not None):
      if (isinstance(parameters, str)):
        parameters = [float(p) for p in parameters.replace(",", " ").split()]
      ext.unit_cell.__init__(self,
        parameters=parameters)
    elif (metrical_matrix is not None):
      ext.unit_cell.__init__(self,
        metrical_matrix=metrical_matrix)
    elif (orthogonalization_matrix is not None):
      ext.unit_cell.__init__(self,
        orthogonalization_matrix=orthogonalization_matrix)
    else:
      ext.unit_cell.__init__(self, parameters=[])

def update_docstring(obj, string):
  try:
    obj.__doc__ = string # Py3
  except AttributeError:
    obj.__func__.__doc__ = string # Py2

update_docstring(ext.unit_cell.metrical_matrix, """
Access to the metrical matrix:

.. math::
   \\begin{pmatrix}
     a^2                & a b \\cos(\\gamma) & a c \\cos(\\beta)  \\\\
     a b \\cos(\\gamma) & b^2                & b c \\cos(\\alpha) \\\\
     a c \\cos(\\beta)  & b c \\cos(\\alpha) & c^2                \\\\
   \\end{pmatrix}

:returns: the metrical matrix
:rtype: tuple
""")

update_docstring(ext.unit_cell.orthogonalization_matrix, """
Matrix for the conversion of fractional to Cartesian coordinates:

.. math::
  \\mathbf{x}_\\textrm{Cartesian} = \\mathbf{O} \\mathbf{x}_\\textrm{fractional}

:returns: the orthogonalization matrix
:rtype: tuple
""")

update_docstring(ext.unit_cell.fractionalization_matrix, """
Matrix for the conversion of fractional to Cartesian coordinates:

.. math::
  \\mathbf{x}_\\textrm{fractional} = \\mathbf{F} \\mathbf{x}_\\textrm{Cartesian}

:returns: the fractionalization matrix
:rtype: tuple
""")

update_docstring(ext.unit_cell.orthogonalize, """
Conversion of fractional coordinates to Cartesian coordinates.

:param sites_frac: The fractional coordinates.
                   Either coordinates for a single site (3-tuple) or a
                   flex.vec3_double array of coordinates.

:returns: the Cartesian coordinates
:rtype: 3-tuple or flex.vec3_double
""")

update_docstring(ext.unit_cell.fractionalize, """
Conversion of Cartesian coordinates to fractional coordinates.

:param sites_cart: The Cartesian coordinates.
                   Either coordinates for a single site (3-tuple) or a
                   flex.vec3_double array of coordinates.

:returns: the fractional coordinates
:rtype: 3-tuple or flex.vec3_double
""")

@bp.inject_into(ext.unit_cell)
class _():

  def __str__(self):
    return format(self, "({:.6g}, {:.6g}, {:.6g}, {:.6g}, {:.6g}, {:.6g})")

  def __repr__(self):
    return format(self, "uctbx.unit_cell(({}, {}, {}, {}, {}, {}))")

  def __format__(self, format_spec="({:.6g}, {:.6g}, {:.6g}, {:.6g}, {:.6g}, {:.6g})"):
    return format_spec.format(*self.parameters())

  def show_parameters(self, f=None, prefix="Unit cell: "):
    if (f is None): f = sys.stdout
    print(prefix + str(self), file=f)

  def is_conventional_hexagonal_basis(self,
        absolute_length_tolerance=1e-3,
        absolute_angle_tolerance=1e-3):
    p = self.parameters()
    if (abs(p[0]-p[1]) > absolute_length_tolerance):
      return False
    ideal_angles = [90, 90, 120]
    for i in [3,4,5]:
      if (abs(p[i]-ideal_angles[i-3]) > absolute_angle_tolerance):
        return False
    return True

  def is_conventional_rhombohedral_basis(self,
        absolute_length_tolerance=1e-3,
        absolute_angle_tolerance=1e-3):
    p = self.parameters()
    for j in [1,2]:
      if (abs(p[0]-p[j]) > absolute_length_tolerance):
        return False
    for j in [4,5]:
      if (abs(p[3]-p[j]) > absolute_angle_tolerance):
        return False
    return True

  def distance_mod_1(self, site_frac_1, site_frac_2):
    return distance_mod_1(
      unit_cell=self, site_frac_1=site_frac_1, site_frac_2=site_frac_2)

  def minimum_reduction(self, iteration_limit=None,
                              multiplier_significant_change_test=None,
                              min_n_no_significant_change=None):
    if (iteration_limit is None):
      iteration_limit = 100
    if (multiplier_significant_change_test is None):
      multiplier_significant_change_test = 16
    if (min_n_no_significant_change is None):
      min_n_no_significant_change = 2
    return fast_minimum_reduction(self,
      iteration_limit,
      multiplier_significant_change_test,
      min_n_no_significant_change)

  def minimum_cell(self, iteration_limit=None,
                         multiplier_significant_change_test=None,
                         min_n_no_significant_change=None):
    return self.minimum_reduction(
      iteration_limit,
      multiplier_significant_change_test,
      min_n_no_significant_change).as_unit_cell()

  def is_buerger_cell(self, relative_epsilon=None):
    from cctbx.uctbx.reduction_base import gruber_parameterization
    return gruber_parameterization(self, relative_epsilon).is_buerger_cell()

  def is_niggli_cell(self, relative_epsilon=None):
    from cctbx.uctbx.reduction_base import gruber_parameterization
    return gruber_parameterization(self, relative_epsilon).is_niggli_cell()

  def niggli_reduction(self, relative_epsilon=None, iteration_limit=None):
    from cctbx.uctbx import krivy_gruber_1976
    return krivy_gruber_1976.reduction(self, relative_epsilon, iteration_limit)

  def niggli_cell(self,
        relative_epsilon=None,
        iteration_limit=None):
    return self.niggli_reduction(
      relative_epsilon, iteration_limit).as_unit_cell()

  def change_of_basis_op_to_niggli_cell(self,
        relative_epsilon=None,
        iteration_limit=None):
    return self.niggli_reduction(
      relative_epsilon, iteration_limit).change_of_basis_op()

  def lattice_symmetry_group(self,
        max_delta=3,
        enforce_max_delta_for_generated_two_folds=True):
    from cctbx import sgtbx
    return sgtbx.lattice_symmetry_group(
      reduced_cell=self,
      max_delta=max_delta,
      enforce_max_delta_for_generated_two_folds
        =enforce_max_delta_for_generated_two_folds)

  def complete_miller_set_with_lattice_symmetry(self,
        anomalous_flag,
        d_min,
        lattice_symmetry_max_delta=3):
    cb_op = self.change_of_basis_op_to_niggli_cell()
    niggli_cell = self.change_basis(cb_op)
    lattice_group = niggli_cell.lattice_symmetry_group(
      max_delta=lattice_symmetry_max_delta)
    from cctbx import crystal
    niggli_lattice_symmetry = crystal.symmetry(
      unit_cell=niggli_cell,
      space_group_info=lattice_group.info())
    from cctbx import miller
    niggli_lattice_set = miller.build_set(
      crystal_symmetry=niggli_lattice_symmetry,
      anomalous_flag=anomalous_flag,
      d_min=d_min)
    return niggli_lattice_set.change_basis(cb_op.inverse())

  def buffer_shifts_frac(self, buffer):
    from cctbx.crystal import direct_space_asu
    return direct_space_asu.float_asu(
      unit_cell=self,
      cuts=[direct_space_asu.float_cut_plane(n=n, c=0)
        for n in [(-1,0,0),(0,-1,0),(0,0,-1)]]) \
      .add_buffer(thickness=float(buffer)) \
      .shape_vertices().max()

  def box_frac_around_sites(self,
        sites_cart=None,
        sites_frac=None,
        buffer=None):
    assert [sites_cart, sites_frac].count(None) == 1
    if (sites_frac is None):
      assert sites_cart.size() > 0
      sites_frac = self.fractionalize(sites_cart=sites_cart)
    else:
      assert sites_frac.size() > 0
    del sites_cart
    if (buffer is None or buffer == 0):
      return sites_frac.min(), sites_frac.max()
    s_min, s_max = sites_frac.min(), sites_frac.max()
    del sites_frac
    shifts_frac = self.buffer_shifts_frac(buffer=buffer)
    return tuple([s-b for s,b in zip(s_min, shifts_frac)]), \
           tuple([s+b for s,b in zip(s_max, shifts_frac)])

  def debye_waller_factors(self,
        miller_index=None,
        miller_indices=None,
        u_iso=None,
        b_iso=None,
        u_cart=None,
        b_cart=None,
        u_cif=None,
        u_star=None,
        exp_arg_limit=50,
        truncate_exp_arg=False):
    assert [miller_index, miller_indices].count(None) == 1
    assert [u_iso, b_iso, u_cart, b_cart, u_cif, u_star].count(None) == 5
    from cctbx import adptbx
    h = miller_index
    if (h is None): h = miller_indices
    if (u_iso is not None):
      b_iso = adptbx.u_as_b(u_iso)
    if (b_iso is not None):
      return adptbx.debye_waller_factor_b_iso(
        self.stol_sq(h),
        b_iso, exp_arg_limit, truncate_exp_arg)
    if (b_cart is not None):
      u_cart = adptbx.b_as_u(b_cart)
    if (u_cart is not None):
      u_star = adptbx.u_cart_as_u_star(self, u_cart)
    if (u_cif is not None):
      u_star = adptbx.u_cif_as_u_star(self, u_cif)
    assert u_star is not None
    return adptbx.debye_waller_factor_u_star(
      h, u_star, exp_arg_limit, truncate_exp_arg)

  debye_waller_factor = debye_waller_factors

def non_crystallographic_buffer_layer(
      sites_cart_min,
      sites_cart_max,
      default_buffer_layer=0.5):
  sites_span = matrix.col(sites_cart_max) - matrix.col(sites_cart_min)
  buffer_layer = max(list(sites_span))
  if (buffer_layer == 0):
    buffer_layer = default_buffer_layer
  return buffer_layer

def non_crystallographic_unit_cell(
      sites_cart=None,
      sites_cart_min=None,
      sites_cart_max=None,
      buffer_layer=None,
      default_buffer_layer=0.5,
      min_unit_cell_length=0):
  assert (sites_cart is None) is not (sites_cart_min is None)
  assert (sites_cart_min is None) is (sites_cart_max is None)
  if (sites_cart is not None):
    sites_cart_min = sites_cart.min()
    sites_cart_max = sites_cart.max()
  if (buffer_layer is None):
    buffer_layer = non_crystallographic_buffer_layer(
      sites_cart_min=sites_cart_min,
      sites_cart_max=sites_cart_max,
      default_buffer_layer=default_buffer_layer)
  sites_span = matrix.col(sites_cart_max) - matrix.col(sites_cart_min)
  return unit_cell([max(min_unit_cell_length, unit_cell_length)
    for unit_cell_length in (sites_span + matrix.col([buffer_layer]*3)*2)])

class non_crystallographic_unit_cell_with_the_sites_in_its_center(object):

  def __init__(self,
        sites_cart,
        buffer_layer=None,
        default_buffer_layer=0.5,
        min_unit_cell_length=0):
    sites_cart_min = sites_cart.min()
    sites_cart_max = sites_cart.max()
    self.unit_cell = non_crystallographic_unit_cell(
      sites_cart=None,
      sites_cart_min=sites_cart_min,
      sites_cart_max=sites_cart_max,
      buffer_layer=buffer_layer,
      default_buffer_layer=default_buffer_layer,
      min_unit_cell_length=min_unit_cell_length)
    unit_cell_center = matrix.col(self.unit_cell.orthogonalize([0.5]*3))
    model_center = (  matrix.col(sites_cart.max())
                    + matrix.col(sites_cart.min())) / 2
    self.shift_vector = unit_cell_center - model_center
    self.sites_cart = sites_cart + self.shift_vector

  def crystal_symmetry(self):
    from cctbx import crystal
    from cctbx import sgtbx
    return crystal.symmetry(
      unit_cell=self.unit_cell,
      space_group=sgtbx.space_group())

  def sites_frac(self):
    return self.unit_cell.fractionalize(self.sites_cart)

def infer_unit_cell_from_symmetry(params, space_group):
  # XXX exercised by iotbx/kriber/tst_strudat.py
  # XXX TODO: add to uctbx tests
  from cctbx import sgtbx
  error_msg = "Cannot interpret unit cell parameters."
  #
  laue_group = str(sgtbx.space_group_info(
    group=space_group.build_derived_laue_group())).replace(" ", "")
  #
  if (len(params) == 6):
    return unit_cell(params)
  else:
    crystal_system = space_group.crystal_system()
    if (crystal_system == "Cubic"):
      if len(params) == 1:
        a = params[0]
      elif len(params) == 3:
        a,b,c = params
        assert a==b==c
      else:
        raise RuntimeError(error_msg)
      unit_cell_ = unit_cell((a,a,a,90,90,90))
    elif (crystal_system in ("Hexagonal", "Trigonal")):
      is_rhombohedral = False
      if (crystal_system == "Trigonal"):
        if (laue_group in ("R-3m:R", "R-3:R")):
          is_rhombohedral = True
      if (is_rhombohedral):
        if len(params) != 2: raise RuntimeError(error_msg)
        a = params[0]
        angle = params[1]
        unit_cell_ = unit_cell((a,a,a,angle,angle,angle))
      else:
        if len(params) == 2:
          a = params[0]
          c = params[1]
        elif len(params) == 3:
          a,b,c = params
          assert a == b
        elif len(params) == 4:
          a,b,c,angle = params
          assert a == b
          assert angle == 120
        else:
          raise RuntimeError(error_msg)
        unit_cell_ = unit_cell((a,a,c,90,90,120))
    elif (crystal_system == "Tetragonal"):
      if len(params) == 2:
        a = params[0]
        c = params[1]
        unit_cell_ = unit_cell((a,a,c,90,90,90))
      elif len(params) == 3:
        a,b,c = params[:3]
        assert a == b
        unit_cell_ = unit_cell((a,a,c,90,90,90))
      else:
        raise RuntimeError(error_msg)
    elif (crystal_system == "Orthorhombic"):
      if len(params) != 3: raise RuntimeError(error_msg)
      a = params[0]
      b = params[1]
      c = params[2]
      unit_cell_ = unit_cell((a,b,c,90,90,90))
    elif (crystal_system == "Monoclinic"):
      if len(params) != 4: raise RuntimeError(error_msg)
      a = params[0]
      b = params[1]
      c = params[2]
      angle = params[3]
      if (laue_group == "P12/m1"):
        unit_cell_ = unit_cell((a,b,c,90,angle,90))
      elif (laue_group == "P112/m"):
        unit_cell_ = unit_cell((a,b,c,90,90,angle))
      elif (laue_group == "P2/m11"):
        unit_cell_ = unit_cell((a,b,c,angle,90,90))
    elif (crystal_system == "Triclinic"):
      raise RuntimeError(error_msg)
    return unit_cell_
