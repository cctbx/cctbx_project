import scitbx.array_family.flex

import boost.python
ext = boost.python.import_ext("scitbx_lbfgs_ext")
from scitbx_lbfgs_ext import *

from scitbx.python_utils.misc import adopt_init_args

class core_parameters:

  def __init__(self, m=5,
                     maxfev=20,
                     gtol=0.9,
                     xtol=1.e-16,
                     stpmin=1.e-20,
                     stpmax=1.e20):
    adopt_init_args(self, locals())

class termination_parameters:

  def __init__(self, traditional_convergence_test=0001,
                     traditional_convergence_test_eps=1.e-5,
                     min_iterations=0,
                     max_iterations=None,
                     max_calls=None):
    adopt_init_args(self, locals())

class exception_handling_parameters:

  def __init__(self, ignore_line_search_failed_rounding_errors=0001,
                     ignore_line_search_failed_step_at_lower_bound=00000,
                     ignore_line_search_failed_step_at_upper_bound=00000,
                     ignore_line_search_failed_maxfev=00000,
                     ignore_line_search_failed_xtol=00000,
                     ignore_search_direction_not_descent=00000):
    adopt_init_args(self, locals())

  def filter(self, msg, n, x, g):
    if (not msg.startswith("lbfgs error")):
      return 1
    if (msg.find("Rounding errors prevent further progress.") >= 0):
      if (self.ignore_line_search_failed_rounding_errors):
        return 0
    elif (msg.find("The step is at the lower bound stpmax().") >= 0):
      if (ext.traditional_convergence_test(n)(x, g)):
        return 0
      if (self.ignore_line_search_failed_step_at_lower_bound):
        return -1
    elif (msg.find("The step is at the upper bound stpmax().") >= 0):
      if (self.ignore_line_search_failed_step_at_upper_bound):
        return -1
    elif (msg.find("Number of function evaluations has reached"
                 + " maxfev().") >= 0):
      if (self.ignore_line_search_failed_maxfev):
        return -1
    elif (msg.find("Relative width of the interval of uncertainty"
                 + " is at most xtol().") >= 0):
      if (self.ignore_line_search_failed_xtol):
        return -1
    elif (msg.find("The search direction is not a descent direction.") >= 0):
      if (ext.traditional_convergence_test(n)(x, g)):
        return 0
      if (self.ignore_search_direction_not_descent):
        return -1
    return 1

def run_c_plus_plus(target_evaluator,
                    termination_params=None,
                    core_params=None,
                    exception_handling_params=None):
  if (termination_params is None):
    termination_params = termination_parameters()
  if (core_params is None):
    core_params = core_parameters()
  if (exception_handling_params is None):
    exception_handling_params = exception_handling_parameters()
  minimizer = ext.minimizer(
    target_evaluator.n,
    core_params.m,
    core_params.maxfev,
    core_params.gtol,
    core_params.xtol,
    core_params.stpmin,
    core_params.stpmax)
  if (termination_params.traditional_convergence_test):
    is_converged = ext.traditional_convergence_test(
      target_evaluator.n,
      termination_params.traditional_convergence_test_eps)
  else:
    is_converged = ext.drop_convergence_test(termination_params.min_iterations)
  callback_after_step = getattr(target_evaluator, "callback_after_step", None)
  x_after_step = None
  x = None
  try:
    while 1:
      x, f, g = target_evaluator()
      if (minimizer.run(x, f, g)): continue
      x_after_step = x.deep_copy()
      if (callback_after_step is not None):
        callback_after_step(minimizer)
      if (termination_params.traditional_convergence_test):
        if (    minimizer.iter() >= termination_params.min_iterations
            and is_converged(x, g)):
          break
      else:
        if (is_converged(f)): break
      if (    termination_params.max_iterations is not None
          and minimizer.iter() >= termination_params.max_iterations):
        break
      if (    termination_params.max_calls is not None
          and minimizer.nfun() > termination_params.max_calls):
        break
      if (not minimizer.run(x, f, g)): break
  except RuntimeError, e:
    if (x is not None and x_after_step is not None):
      x.clear()
      x.extend(x_after_step)
    minimizer.error = str(e)
    error_classification = exception_handling_params.filter(
      minimizer.error, target_evaluator.n, x, g)
    if (error_classification > 0):
      raise
    elif (error_classification < 0):
      minimizer.is_unusual_error = 0001
    else:
      minimizer.is_unusual_error = 00000
  else:
    minimizer.error = None
    minimizer.is_unusual_error = None
  return minimizer

def run_fortran(target_evaluator,
                termination_params=None,
                core_params=None):
  "For debugging only!"
  from scitbx.python_utils.misc import store
  from fortran_lbfgs import lbfgs as fortran_lbfgs
  import Numeric
  if (termination_params is None):
    termination_params = termination_parameters()
  if (core_params is None):
    core_params = core_parameters()
  assert termination_params.traditional_convergence_test
  assert core_params.maxfev == 20
  n = target_evaluator.n
  m = core_params.m
  gtol = core_params.gtol
  xtol = core_params.xtol
  stpmin = core_params.stpmin
  stpmax = core_params.stpmax
  eps = termination_params.traditional_convergence_test_eps
  x = Numeric.array(Numeric.arange(n), Numeric.Float64)
  g = Numeric.array(Numeric.arange(n), Numeric.Float64)
  size_w = n*(2*m+1)+2*m
  w = Numeric.array(Numeric.arange(size_w), Numeric.Float64)
  diag = Numeric.array(Numeric.arange(n), Numeric.Float64)
  iprint = [1, 0]
  diagco = 0
  iflag = Numeric.array([0], Numeric.Int32)
  minimizer = store(error=None)
  while 1:
    x_, f, g_ = target_evaluator()
    for i,xi in enumerate(x_): x[i] = xi
    for i,gi in enumerate(g_): g[i] = gi
    fortran_lbfgs(n, m, x, f, g, diagco, diag, iprint, eps, xtol, w, iflag)
    for i,xi in enumerate(x): x_[i] = xi
    if (iflag[0] == 0):
      break
    if (iflag[0] < 0):
      minimizer.error = "fortran lbfgs error"
      break
  return minimizer

def run(target_evaluator,
        termination_params=None,
        core_params=None,
        exception_handling_params=None,
        use_fortran=00000):
  if (use_fortran):
    return run_fortran(target_evaluator, termination_params, core_params)
  else:
    return run_c_plus_plus(
      target_evaluator,
      termination_params,
      core_params,
      exception_handling_params)
