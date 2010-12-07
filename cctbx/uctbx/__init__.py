import cctbx.array_family.flex # import dependency
import scitbx.array_family.shared # import dependency

import boost.python
ext = boost.python.import_ext("cctbx_uctbx_ext")
from cctbx_uctbx_ext import *

from scitbx import matrix
import sys

class unit_cell(ext.unit_cell):

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

class _unit_cell(boost.python.injector, ext.unit_cell):

  def __str__(self):
    return "(%.6g, %.6g, %.6g, %.6g, %.6g, %.6g)" % self.parameters()

  def show_parameters(self, f=None, prefix="Unit cell: "):
    if (f is None): f = sys.stdout
    print >> f, prefix + str(self)

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

  def niggli_cell(self, relative_epsilon=None, iteration_limit=None):
    return self.niggli_reduction(
      relative_epsilon, iteration_limit).as_unit_cell()

  def buffer_shifts_frac(self, buffer):
    from cctbx.crystal import direct_space_asu
    return direct_space_asu.float_asu(
      unit_cell=self,
      cuts=[direct_space_asu.float_cut_plane(n=n, c=0)
        for n in [(-1,0,0),(0,-1,0),(0,0,-1)]]) \
      .add_buffer(thickness=float(buffer)) \
      .volume_vertices().max()

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
    self.sites_cart = sites_cart + (unit_cell_center - model_center)

  def crystal_symmetry(self):
    from cctbx import crystal
    from cctbx import sgtbx
    return crystal.symmetry(
      unit_cell=self.unit_cell,
      space_group=sgtbx.space_group())
