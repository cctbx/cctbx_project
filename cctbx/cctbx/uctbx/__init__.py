from cctbx.array_family import flex
import scitbx.array_family.shared

import boost.python
ext = boost.python.import_ext("cctbx_uctbx_ext")
from cctbx_uctbx_ext import *

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

  def minimum_reduction(self, iteration_limit=None,
                              multiplier_significant_change_test=None,
                              min_n_no_significant_change=None):
    if (iteration_limit is None):
      iteration_limit = 100
    if (multiplier_significant_change_test is None):
      multiplier_significant_change_test = 10
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

  def box_frac_around_sites(self,
        sites_cart=None,
        sites_frac=None,
        buffer=None):
    assert [sites_cart, sites_frac].count(None) == 1
    if (buffer is None):
      if (sites_frac is None):
        assert sites_cart.size() > 0
        sites_frac = self.fractionalization_matrix() * sites_cart
      else:
        assert sites_frac.size() > 0
      del sites_cart
      return sites_frac.min(), sites_frac.max()
    assert isinstance(buffer, float)
    if (sites_cart is None):
      assert sites_frac.size() > 0
      sites_cart = self.orthogonalization_matrix() * sites_frac
    else:
      assert sites_cart.size() > 0
    del sites_frac
    min_max = []
    for b in [-buffer, buffer]:
      sites_frac = self.fractionalization_matrix() * (sites_cart + [b]*3)
      min_max.append([sites_frac.min(), sites_frac.max()])
    return tuple([min(x,y) for x,y in zip(min_max[0][0], min_max[1][0])]), \
           tuple([max(x,y) for x,y in zip(min_max[0][1], min_max[1][1])])
