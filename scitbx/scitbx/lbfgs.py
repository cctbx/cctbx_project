import scitbx.array_family.flex

from scitbx.python_utils import misc
ext = misc.import_ext("scitbx_boost.lbfgs_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

def run(target_evaluator,
        min_iterations=10,
        max_calls=100,
        traditional_convergence_test=0):
  minimizer = ext.minimizer(target_evaluator.n)
  if (traditional_convergence_test):
    is_converged = ext.traditional_convergence_test(target_evaluator.n)
  else:
    is_converged = ext.drop_convergence_test(min_iterations)
  try:
    while 1:
      x, f, g = target_evaluator()
      if (minimizer.run(x, f, g)): continue
      if (traditional_convergence_test):
        if (minimizer.iter() >= min_iterations and is_converged(x, g)): break
      else:
        if (is_converged(f)): break
      if (minimizer.nfun() > max_calls): break
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
  return minimizer
