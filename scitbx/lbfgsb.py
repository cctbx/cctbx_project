from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("scitbx_lbfgsb_ext")
from scitbx_lbfgsb_ext import *


class minimizer(ext.minimizer):

  def __init__(self, n=None,
                     m=None,
                     l=None,
                     u=None,
                     nbd=None,
                     enable_stp_init=False,
                     factr=None,
                     pgtol=None,
                     iprint=None):
    assert [l,u,nbd].count(None) in [0,3]
    assert n is not None or l is not None
    if (n is None):
      assert u.size() == l.size() and nbd.size() == l.size()
      n = l.size()
    elif (l is None):
      l = flex.double(n, 0)
      u = l
      nbd = flex.int(n, 0)
    if (m is None): m = 5
    if (factr is None): factr = 1.0e+7
    if (pgtol is None): pgtol = 1.0e-5
    if (iprint is None): iprint = -1
    ext.minimizer.__init__(self,
      n, m, l, u, nbd, enable_stp_init, factr, pgtol, iprint)


def run(target_evaluator,
               bound_flags,
               lower_bound,
               upper_bound,
               n,
               max_iterations = None):
  """
    Run the L-BFGS-B minimization algorithm using the provided target
    evaluator and bounds.

    This function minimizes a target function using the Limited-memory
    Broyden-Fletcher-Goldfarb-Shanno with Box constraints (L-BFGS-B)
    algorithm. It iteratively calls the target evaluator to compute the
    function value and gradients, and processes them with the L-BFGS-B
    minimizer.

    Parameters
    ----------
    target_evaluator : callable
        A callable object that provides the method
        `compute_functional_and_gradients()`, which returns the current
        values of the variables `x`, the function value `f`, and the
        gradient `g`.
    use_bounds : int
        Number of bounds used in the optimization. Used to compute nbd which
        is the flag that tells the optimizer which bounds to use.
    lower_bound : flex.double
        The lower bound for the optimization variables.
    upper_bound : flex.double
        The upper bound for the optimization variables.
    n : int
        The size of x (i.e. the number of variables to be optimized).
    max_iterations : int, optional
        Maximum number of iterations allowed during the optimization. If
        None, there is no limit on the number of iterations.

    Returns
    -------
    lbfgsb_minimizer : object
        An instance of the minimizer object that contains the result of the
        optimization process. Its attributes are:
        - `x`: The optimized variables.
        - `f`: The minimized function value.
        - `g`: The gradient at the optimized point.
        - `error`: A string containing any error message that occurred
          during execution.
        - `n_calls`: The number of iterations performed.

    Raises
    ------
    RuntimeError
        If the L-BFGS-B minimizer encounters an error during execution. The
        error message is stored in `lbfgsb_minimizer.error`.

    Notes
    -----
    - The process runs in a loop until either the minimizer converges, the
      specified maximum number of iterations is reached, or an error
      occurs.
    - If bounds are specified, the optimization is constrained within the
      provided lower and upper bounds.
  """
  callback_after_step = getattr(target_evaluator, "callback_after_step", None)
  lbfgsb_minimizer = minimizer(
    n   = n,
    l   = lower_bound,
    u   = upper_bound,
    nbd = bound_flags) # flag to apply both bounds
  lbfgsb_minimizer.error = None
  x, f, g = target_evaluator.compute_functional_and_gradients()
  try:
    icall = 0
    while 1:
      have_request = lbfgsb_minimizer.process(x, f, g)
      if(have_request):
        requests_f_and_g = lbfgsb_minimizer.requests_f_and_g()
        x, f, g = target_evaluator.compute_functional_and_gradients()
        icall += 1
        continue
      assert not lbfgsb_minimizer.requests_f_and_g()
      if(lbfgsb_minimizer.is_terminated()): break
      if(max_iterations is not None and icall>max_iterations): break
      if(callback_after_step is not None):
        if(callback_after_step(minimizer) is True):
          print("lbfgs minimizer stop: callback_after_step is True")
          break
  except RuntimeError as e:
    lbfgsb_minimizer.error = str(e)
  lbfgsb_minimizer.n_calls = icall
  if(lbfgsb_minimizer.error is not None):
    print("lbfgs-b: an error occured: %s"%lbfgsb_minimizer.error)
  return lbfgsb_minimizer
