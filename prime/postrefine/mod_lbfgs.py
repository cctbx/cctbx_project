from __future__ import absolute_import, division, print_function
from scitbx import lbfgs
from cctbx.array_family import flex
from .mod_lbfgs_partiality import lbfgs_partiality_handler
from six.moves import range

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
    lp_h = lbfgs_partiality_handler()
    #calculate sum_sqr of the function
    fvec = lp_h.func(self.x, self.args)
    self.f = flex.sum(fvec*fvec)
    #calculate gradient for each parameter
    DELTA = 1.E-7
    self.g = flex.double()
    for x in range(self.n):
      templist = list(self.x)
      templist[x]+=DELTA
      dvalues = flex.double(templist)
      dfvec = lp_h.func(dvalues, self.args)
      df = flex.sum(dfvec*dfvec)
      #calculate by finite_difference
      self.g.append( ( df-self.f )/DELTA )
    return self.f, self.g
