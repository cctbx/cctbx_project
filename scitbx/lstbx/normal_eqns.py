import boost.python
boost.python.import_ext("scitbx_lstbx_normal_equations_ext")
import scitbx_lstbx_normal_equations_ext as ext
from scitbx_lstbx_normal_equations_ext import *

class non_linear_ls_mixin(object):
  """ Synopsis:

  Either of

  class non_linear_ls_xxx(core_non_linear_ls_xxx,
                          non_linear_ls_mixin):
    pass

  class non_linear_ls_xxx(non_linear_ls_mixin):

    # define required methods

  This class is a mixin to be inherited as shown. It requires
  core_non_linear_ls_xxx to provide some methods in the first case,
  or those methods to be defined directly in the heir in the second case:
  those that this mixin marks as not-implemented.
  """

  def parameter_vector_norm(self):
    raise NotImplementedError()

  def build_up(self, objective_only=False):
    raise NotImplementedError()

  def step_equations(self):
    raise NotImplementedError()

  def objective(self):
    raise NotImplementedError()

  def step_forward(self):
    raise NotImplementedError()

  def opposite_of_gradient(self):
    return self.step_equations().right_hand_side()

  def normal_matrix_packed_u(self):
    return self.step_equations().normal_matrix_packed_u()

  def solve(self):
    self.step_equations().solve()

  def solve_and_step_forward(self):
    self.solve()
    self.step_forward()

  def step(self):
    return self.step_equations().solution()


class _(boost.python.injector, linear_ls):

  def __iter__(self):
    yield self.normal_matrix_packed_u()
    yield self.right_hand_side()


class non_linear_ls_with_separable_scale_factor(
  ext.non_linear_ls_with_separable_scale_factor,
  non_linear_ls_mixin):

  pass
