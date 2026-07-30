from __future__ import absolute_import, division, print_function
import libtbx.load_env
import boost_adaptbx.boost.python as bp
bp.import_ext("scitbx_lstbx_normal_equations_ext")
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

  def apply_shifts(self, shifts):
    """ Move the parameter vector by the given increment.

    step_forward() is the special case shifts = step(). Minimisers which do not
    follow the step obtained from the normal equations -- e.g. those in
    lstbx.scipy_iterations -- need to move to a point of their own choosing and
    therefore require this instead.
    """
    raise NotImplementedError()

  def parameter_vector(self):
    """ The current parameter vector, or None if the problem cannot report it.

    Only needed by minimisers which move the parameters to arbitrary points:
    knowing where the parameters actually ended up lets them recover from
    values that were clamped while the shifts were applied. Returning None is
    always safe; such a minimiser then keeps track of the shifts it applied.
    """
    return None

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

  def save_step_equations(self):
    """ Copies of the step equations, to be put back with restore later.

    Solving factorises the normal matrix in place and overwrites the right hand
    side with the solution, so equations cannot be solved twice -- with a
    different Levenberg-Marquardt damping, say. Rebuilding them costs a pass
    over the data; copying them beforehand costs the matrix, once, and no data
    are touched at all.

    Call before the damping is added, so that what comes back is the undamped
    normal equations.
    """
    return (self.normal_matrix_packed_u().deep_copy(),
            self.opposite_of_gradient().deep_copy())

  def restore_step_equations(self, saved):
    """ Put back what save_step_equations returned, unsolved.

    reset() zeroes both arrays and clears the solved flag; the arrays it hands
    out share the storage the equations themselves use, so adding the copies
    back in leaves the object as it was before it was damped and solved.
    """
    normal_matrix, right_hand_side = saved
    equations = self.step_equations()
    equations.reset()
    # bound to a name first: the arrays share the equations' own storage, so
    # adding into them in place is what puts the values back
    a = equations.normal_matrix_packed_u()
    a += normal_matrix
    b = equations.right_hand_side()
    b += right_hand_side

@bp.inject_into(linear_ls)
class _():

  def __iter__(self):
    yield self.normal_matrix_packed_u()
    yield self.right_hand_side()


non_linear_ls_with_separable_scale_factor_description = """\
* non-linear L.S. by optimising the overall scale factor alone first
  and then the other parameters alone
"""

class non_linear_ls_with_separable_scale_factor_BLAS_2(
  ext.non_linear_ls_with_separable_scale_factor__level_2_blas_impl,
  non_linear_ls_mixin):

  @property
  def description(self):
    return (non_linear_ls_with_separable_scale_factor_description +
            "* slow normal matrix computation\n")

  @property
  def debug_info(self):
    return ""


if libtbx.env.has_module('fast_linalg'):
  class non_linear_ls_with_separable_scale_factor_BLAS_3(
    ext.non_linear_ls_with_separable_scale_factor__level_3_blas_impl,
    non_linear_ls_mixin):

    @property
    def description(self):
      return (non_linear_ls_with_separable_scale_factor_description +
              "* fast normal matrix computation")

    @property
    def debug_info(self):
      import fast_linalg
      e = fast_linalg.env
      return '\n'.join((
        "\n", "*"*80,
        "*** Using OpenBLAS with %i threads on a machine with %i %s cores" %
        (e.threads, e.physical_cores, e.cpu_family),
        "*** OpenBLAS was built with the following options:",
        "*** %s" % e.build_config,
        "*"*80, "\n",
      ))


non_linear_ls_with_separable_scale_factor = \
  non_linear_ls_with_separable_scale_factor_BLAS_2
