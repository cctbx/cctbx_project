"""
lbfgs.tst_func_free_line_search

"""

from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import scitbx.lbfgs
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal

def func(x1,x2,x3):
  return (x1-2)**2 + (x2+5)**2 + (x3-500)**4

def grad(x1,x2,x3):
  return flex.double([(x1-2)*2, (x2+5)*2, 4*(x3-500)*3])

x0 = [-100,203,999]

class lbfgs(object):

  def __init__(self, gradient_only):
    print("scitbx.lbfgs: use gradient_only line search: ", gradient_only)
    adopt_init_args(self, locals())
    self.lbfgs_core_params = scitbx.lbfgs.core_parameters(
                               stpmin=1.e-10,
                               stpmax=500)
    self.lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = 6000)
    self.f=None
    self.g=None
    self.x = flex.double(x0)
    self.minimizer = scitbx.lbfgs.run(
      gradient_only             = gradient_only,
      target_evaluator          = self,
      core_params               = self.lbfgs_core_params,
      termination_params        = self.lbfgs_termination_params,
      exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_step_at_lower_bound = True))

  def compute_functional_and_gradients(self):
    self.f = func(x1=self.x[0], x2 = self.x[1], x3 = self.x[2])
    self.g = grad(x1=self.x[0], x2 = self.x[1], x3 = self.x[2])
    return self.f, self.g

if(__name__ == "__main__"):
  print("The answer is: [2, -5, 500] ")
  print("The start value is:", x0)
  t=lbfgs(True)
  assert approx_equal(t.x, [2, -5, 500])
  t=lbfgs(False)
  assert approx_equal(t.x, [2, -5, 500])
  print("OK")
