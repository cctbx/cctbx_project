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
                     min_iterations=0,
                     max_iterations=None,
                     max_calls=None):
    adopt_init_args(self, locals())

def run(target_evaluator,
        termination_params=None,
        core_params=None):
  if (termination_params is None):
    termination_params = termination_parameters()
  if (core_params is None):
    core_params = core_parameters()
  minimizer = ext.minimizer(
    target_evaluator.n,
    core_params.m,
    core_params.maxfev,
    core_params.gtol,
    core_params.xtol,
    core_params.stpmin,
    core_params.stpmax)
  if (termination_params.traditional_convergence_test):
    is_converged = ext.traditional_convergence_test(target_evaluator.n)
  else:
    is_converged = ext.drop_convergence_test(termination_params.min_iterations)
  callback_after_step = getattr(target_evaluator, "callback_after_step", None)
  try:
    while 1:
      x, f, g = target_evaluator()
      if (minimizer.run(x, f, g)): continue
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
    minimizer.error = str(e)
    if (    str(e).find("Tolerances may be too small.")<0
        and str(e).find("The search direction is not a descent direction")<0
        and str(e).find("The step is at the lower bound")<0):
      raise
    if (str(e).find("The step is at the lower bound")>=0):
      if (not ext.traditional_convergence_test(target_evaluator.n)(x, g)):
        raise
  else:
    minimizer.error = None
  return minimizer
