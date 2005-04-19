from cctbx.array_family import flex # for tuple mappings

import boost.python
ext = boost.python.import_ext("cctbx_adptbx_ext")
from cctbx_adptbx_ext import *

import scitbx.math
import random

def random_rotate_ellipsoid(u_cart):
  c = scitbx.math.euler_angles_as_matrix(
    [random.uniform(0,360) for i in xrange(3)]).elems
  return c_u_c_transpose(c, u_cart)

def random_u_cart(u_scale=1, u_min=0):
  return random_rotate_ellipsoid(u_cart=[random.random()*u_scale+u_min
    for i in xrange(3)] + [0,0,0])

class constraints:

  def __init__(self, space_group, initialize_gradient_handling=False):
    self.row_echelon_form = space_group.tensor_constraints(
      reciprocal_space=True)
    independent_flags = flex.bool(6, True)
    scitbx.math.row_echelon_back_substitution_int(
      row_echelon_form=self.row_echelon_form,
      independent_flags=independent_flags)
    self.independent_indices = independent_flags.iselection()
    if (not initialize_gradient_handling):
      self.gradient_average_denominator = 0
    else:
      self.gradient_average_denominator = space_group.n_smx()
      self.gradient_average_cache = \
        scitbx.math.tensor_rank_2_gradient_average_cache_int()
      for i in xrange(self.gradient_average_denominator):
        r = space_group(i).r()
        assert r.den() == 1
        self.gradient_average_cache.accumulate(a=r.num())
      self.gradient_sum_coeffs = []
      for indep_index in self.independent_indices:
        coeffs = flex.double(6, 0)
        coeffs[indep_index] = 1
        scitbx.math.row_echelon_back_substitution_float(
          row_echelon_form=self.row_echelon_form,
          solution=coeffs)
        self.gradient_sum_coeffs.append(coeffs)

  def n_independent_params(self):
    return self.independent_indices.size()

  def n_dependent_params(self):
    return 6-self.independent_indices.size()

  def independent_params(self, u_star):
    return flex.double(u_star).select(self.independent_indices)

  def all_params(self, independent_params):
    result = flex.double(6, 0)
    result.set_selected(self.independent_indices, independent_params)
    scitbx.math.row_echelon_back_substitution_float(
      row_echelon_form=self.row_echelon_form,
      solution=result)
    return result

  def sym_gradients(self, asu_gradients):
    if (self.gradient_average_denominator == 0):
      raise RuntimeError("gradient handling was not initialized.")
    return self.gradient_average_cache.average(
      g=asu_gradients,
      denominator=self.gradient_average_denominator)

  def independent_gradients(self, all_gradients, are_sym=False):
    if (not are_sym):
      all_gradients = self.sym_gradients(asu_gradients=all_gradients)
    result = []
    for coeffs in self.gradient_sum_coeffs:
      result.append(flex.sum(coeffs * flex.double(all_gradients)))
    return result
