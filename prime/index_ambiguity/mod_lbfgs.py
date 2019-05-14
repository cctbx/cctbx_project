from __future__ import absolute_import, division, print_function
from scitbx import lbfgs
from cctbx.array_family import flex
from prime import index_ambiguity
import numpy as np

class lbfgs_handler(object):
  """
  lbfgs_handler optimize set of parameters (params) by fitting
  data[0] to data [1] using given function (func).
  """
  def __init__(self, current_x=None, args=None,
               min_iterations=0, max_iterations=None, max_calls=1000, max_drop_eps=1.e-5):
    self.n = current_x.size()
    self.x = current_x
    self.args = args
    self.minimizer = lbfgs.run(
            target_evaluator=self,
            termination_params=lbfgs.termination_parameters(
              traditional_convergence_test=False,
              drop_convergence_test_max_drop_eps=max_drop_eps,
              min_iterations=min_iterations,
              max_iterations=max_iterations,
              max_calls=max_calls),
            exception_handling_params=lbfgs.exception_handling_parameters(
               ignore_line_search_failed_rounding_errors=True,
               ignore_line_search_failed_step_at_lower_bound=True,
               ignore_line_search_failed_step_at_upper_bound=False,
               ignore_line_search_failed_maxfev=False,
               ignore_line_search_failed_xtol=False,
               ignore_search_direction_not_descent=False)
            )

  def compute_functional_and_gradients(self):
    x_set = np.array(self.x).reshape((len(self.x)//2,2))
    r_grid = flex.double(self.args)
    x_vec_set = flex.vec2_double(x_set)
    self.f = index_ambiguity.calc_BD_alg_2_sum_sqr(r_grid, x_vec_set)
    self.g = index_ambiguity.calc_gradient_BD_alg_2(r_grid, x_vec_set)
    return self.f, self.g
