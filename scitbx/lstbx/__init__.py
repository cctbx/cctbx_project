import boost.python
boost.python.import_ext("scitbx_lstbx_ext")
import scitbx_lstbx_ext as ext
from scitbx_lstbx_ext import *

class non_linear_normal_equations_mixin(object):
  """ Synopsis:

  Either of

  class normal_equations_xxx(core_normal_equations_xxx,
                             non_linear_normal_equations_mixin):
    pass

  class normal_equations_xxx(non_linear_normal_equations_mixin):

    # define required methods

  This class is a mixin to be inherited as shown. It requires
  core_normal_equations_xxx to provide some methods in the first case,
  or those methods to be defined directly in the heir in the second case:
  those that this mixin marks as not-implemented.
  """

  def parameter_vector_norm(self):
    raise NotImplementedError()

  def build_up(self, objective_only):
    raise NotImplementedError()

  def step_equations(self):
    raise NotImplementedError()

  def objective(self):
    raise NotImplementedError()

  def step_forward(self):
    raise NotImplementedError()

  def gradient(self):
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


class normal_equations_extension(boost.python.injector, normal_equations):

  def __iter__(self):
    yield self.normal_matrix_packed_u()
    yield self.right_hand_side()


class normal_equations_separating_scale_factor(
  ext.normal_equations_separating_scale_factor,
  non_linear_normal_equations_mixin):

  pass
