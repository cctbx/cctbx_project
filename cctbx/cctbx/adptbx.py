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

class constraints:

  def __init__(self, space_group):
    self.row_echelon_form = space_group.tensor_constraints(
      reciprocal_space=True)
    self.independent_flags = flex.bool(6, True)
    scitbx.math.row_echelon_back_substitution_int(
      row_echelon_form=self.row_echelon_form,
      independent_flags=self.independent_flags)

  def n_independent_params(self):
    return self.independent_flags.count(True)

  def n_dependent_params(self):
    return self.independent_flags.count(False)

  def independent_params(self, u_star):
    return flex.double(u_star).select(self.independent_flags)

  def all_params(self, independent_params):
    result = flex.double(6, 0)
    result.set_selected(self.independent_flags, independent_params)
    scitbx.math.row_echelon_back_substitution_float(
      row_echelon_form=self.row_echelon_form,
      solution=result)
    return result
