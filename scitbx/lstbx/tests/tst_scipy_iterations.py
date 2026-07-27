from __future__ import absolute_import, division, print_function

from scitbx.array_family import flex
from scitbx.lstbx import normal_eqns, normal_eqns_solving
from scitbx.lstbx.tests import test_problems
from libtbx.test_utils import approx_equal
import math

non_linear_ls_with_separable_scale_factor__impls = (
  normal_eqns.non_linear_ls_with_separable_scale_factor_BLAS_2,
)

# Tolerances tight enough that a method stops on its own rather than on the
# iteration cap. Asking for much more than this is asking for more than double
# precision affords on these problems, and a minimiser told to keep going past
# the point where its own gradient has gone to zero does not always survive it.
converged = {'n_max_iterations': 2000, 'max_evaluations': 50000,
             'g_tolerance': 1e-10, 'f_tolerance': 1e-15, 'x_tolerance': 1e-10}


def objective_and_gradient_of(problem):
  """ The objective and its gradient as plain functions of x.

  Only usable on the test problems, whose parameter vector is a bare attribute.
  """
  def objective(x):
    problem.x = flex.double(x)
    problem.build_up()
    return problem.objective()
  def gradient(x):
    problem.x = flex.double(x)
    problem.build_up()
    return (-problem.opposite_of_gradient()).as_numpy_array()
  return objective, gradient


def exercise_gradient_consistency():
  """ opposite_of_gradient() shall really be minus the gradient.

  Everything else here rests on this. The solvers of normal_eqns_solving never
  do a line search and are therefore forgiving of a gradient which is slightly
  off, whereas the scipy ones stall on one. The separable scale factor case is
  the one worth pinning down, its gradient carrying a term through dK*/dx.
  """
  from scipy import optimize
  for impl in non_linear_ls_with_separable_scale_factor__impls:
    for problem, x, name in (
      (test_problems.polynomial_fit(impl)(normalised=True),
       flex.double((0.5, 0.3, 0.2)), 'polynomial_fit'),
      (test_problems.polynomial_fit_with_penalty(impl)(normalised=True),
       flex.double((0.5, 0.3, 0.2)), 'polynomial_fit_with_penalty'),
      (test_problems.exponential_fit(),
       flex.double((-1, -2, 1, -1)), 'exponential_fit')):
      objective, gradient = objective_and_gradient_of(problem)
      error = optimize.check_grad(objective, gradient, x.as_numpy_array())
      norm = flex.double(gradient(x.as_numpy_array())).norm()
      assert error <= 1e-5*max(norm, 1), (name, error, norm)


def exercise_hessian_product():
  """ The Hessian product shall be the normal matrix in preconditioned units.

  The Newton-type methods are handed products with the Gauss-Newton Hessian
  rather than the matrix itself, so nothing else would notice were the
  preconditioner applied on one side only, or the wrong way round.
  """
  import numpy
  from scitbx.lstbx import scipy_iterations
  for impl in non_linear_ls_with_separable_scale_factor__impls:
    problem = test_problems.polynomial_fit(impl)(normalised=True)
    # n_max_iterations=0 builds the normal equations and stops there
    cycles = scipy_iterations.scipy_iterations(
      problem, method='Newton-CG', n_max_iterations=0)
    n = problem.x.size()
    a = cycles.non_linear_ls.normal_matrix_packed_u()
    a = a.matrix_packed_u_as_symmetric().as_numpy_array().reshape(n, n)
    d = numpy.diag(1./cycles.parameter_scale.as_numpy_array())
    expected = d.dot(a).dot(d)
    got = numpy.array([cycles.hessian_product(numpy.zeros(n), e)
                       for e in numpy.eye(n)]).T
    assert approx_equal(got, expected, eps=1e-12), (got, expected)
    # that is what the preconditioner is for
    assert approx_equal(numpy.diag(expected), numpy.ones(n), eps=1e-12)


def exercise_minimisers():
  """ Every supported method shall find the minimum of polynomial_fit.

  Levenberg-Marquardt provides the reference: the data are noisy, so where the
  minimum lies exactly is not known beforehand and agreeing with the solver
  already trusted for this problem is more to the point than agreeing with the
  noiseless answer. The scale factor is optimised away rather than refined, so
  it also has to come out right without ever having been a parameter.
  """
  from scitbx.lstbx import scipy_iterations
  for impl in non_linear_ls_with_separable_scale_factor__impls:
    reference = test_problems.polynomial_fit(impl)(normalised=True)
    normal_eqns_solving.levenberg_marquardt_iterations(
      reference, gradient_threshold=1e-12, step_threshold=1e-12,
      tau=1e-4, n_max_iterations=200)
    assert approx_equal(reference.x, reference.arg_min, eps=1e-3), \
      list(reference.x)

    for method in sorted(scipy_iterations.scipy_iterations.supported_methods):
      problem = test_problems.polynomial_fit(impl)(normalised=True)
      try:
        cycles = scipy_iterations.scipy_iterations(
          problem, method=method, **converged)
      except ImportError as e:
        # trust-krylov needs an optional scipy extension
        print("\tskipped %s (%s)" % (method, e))
        continue
      assert approx_equal(problem.x, reference.x, eps=1e-5), (
        method, list(problem.x), list(reference.x))
      assert approx_equal(problem.optimal_scale_factor(), 2, eps=1e-3), (
        method, problem.optimal_scale_factor())
      print("\t%-13s %4d iterations, %5d evaluations"
            % (method, cycles.n_iterations, cycles.n_evaluations))


def exercise_ill_conditioned_problem():
  """ exponential_fit is hard, and the methods shall differ honestly on it.

  Its minimum sits at the bottom of a very flat valley, so how close a method
  gets to arg_min says more about the method than about whether it is wired up
  correctly. What every method must do is lower the objective; the quasi-Newton
  and trust-region ones are additionally expected to match Levenberg-Marquardt.

  Newton-CG is not among them: its line search works from a truncated CG
  direction, and on a Gauss-Newton Hessian this ill conditioned that direction
  is poor. trust-ncg uses the same Hessian products and does fine, the trust
  region being what makes the difference.
  """
  from scitbx.lstbx import scipy_iterations
  reference = test_problems.exponential_fit()
  normal_eqns_solving.levenberg_marquardt_iterations(
    reference, gradient_threshold=1e-12, step_threshold=1e-12, tau=1e-4,
    n_max_iterations=500)
  assert approx_equal(reference.x, reference.arg_min, eps=5e-4)
  best = reference.objective()

  as_good_as_lm = ('BFGS', 'L-BFGS-B', 'SLSQP', 'trust-ncg')
  for method in sorted(scipy_iterations.scipy_iterations.supported_methods):
    problem = test_problems.exponential_fit()
    problem.build_up()
    start = problem.objective()
    try:
      cycles = scipy_iterations.scipy_iterations(
        problem, method=method, **converged)
    except ImportError:
      continue
    achieved = problem.objective()
    assert achieved < start, (method, achieved, start)
    if method in as_good_as_lm:
      assert achieved <= best*(1 + 1e-6), (method, achieved, best)
      assert approx_equal(problem.x, problem.arg_min, eps=5e-4), (
        method, list(problem.x))
    print("\t%-13s %4d iterations, %5d evaluations, objective %.8e"
          % (method, cycles.n_iterations, cycles.n_evaluations, achieved))


def exercise_preconditioning():
  """ Badly scaled parameters shall be handled, and by the preconditioner.

  Both halves matter: that the minimisers cope, and that it is the scaling of
  the parameters which lets them. Without it the first-order methods run out of
  iterations nowhere near the minimum. trust-ncg is much less sensitive to the
  scaling, a trust region being a length scale of its own, so it is checked to
  converge either way.
  """
  from scitbx.lstbx import scipy_iterations
  for impl in non_linear_ls_with_separable_scale_factor__impls:
    problem = test_problems.badly_scaled_basis_fit(impl)()
    problem.build_up()
    diagonal = problem.normal_matrix_packed_u().matrix_packed_u_diagonal()
    assert flex.max(diagonal)/flex.min(diagonal) > 1e7, \
      flex.max(diagonal)/flex.min(diagonal)

    for method, scaling_needed in (('CG', True), ('L-BFGS-B', True),
                                   ('trust-ncg', False)):
      achieved = {}
      for preconditioning in ('diagonal', None):
        problem = test_problems.badly_scaled_basis_fit(impl)()
        scipy_iterations.scipy_iterations(
          problem, method=method, preconditioning=preconditioning, **converged)
        achieved[preconditioning] = problem.objective()
      assert achieved['diagonal'] < 1e-10, (method, achieved)
      if scaling_needed:
        assert achieved[None] > 1e6*achieved['diagonal'], (method, achieved)
      else:
        assert achieved[None] < 1e-10, (method, achieved)
      print("\t%-10s objective %.3e preconditioned, %.3e without"
            % (method, achieved['diagonal'], achieved[None]))


def exercise_interrupt():
  """ A minimisation abandoned part way shall leave usable parameters. """
  from scitbx.lstbx import scipy_iterations

  class interrupting_iterations(scipy_iterations.scipy_iterations):
    stop_after = 3
    parameters_when_stopped = None

    def on_iteration_completion(self):
      if self.n_iterations >= self.stop_after:
        self.parameters_when_stopped = self.non_linear_ls.x.deep_copy()
        self.interrupted = True

  problem = test_problems.exponential_fit()
  cycles = interrupting_iterations(problem, method='L-BFGS-B',
                                   n_max_iterations=500)
  assert cycles.n_iterations == cycles.stop_after, cycles.n_iterations
  # The parameters are those of the last iterate the minimiser accepted, not of
  # whichever trial point it happened to be evaluating when it was stopped.
  assert approx_equal(problem.x, cycles.parameters_when_stopped, eps=1e-12), (
    list(problem.x), list(cycles.parameters_when_stopped))


def exercise_non_finite_guard():
  """ A minimisation which blows up shall leave finite parameters behind.

  A model overflowing, or a minimiser dividing by a gradient it has itself
  driven to zero, both end in NaN. The parameters are moved by increments, so
  there is no shifting back from one: the guard has to stop them ever going
  there, rather than tidy up afterwards.
  """
  import numpy
  from scitbx.lstbx import scipy_iterations

  class blowing_up(scipy_iterations.scipy_iterations):
    """ Hands the minimiser a NaN gradient once it has taken a few steps. """
    blow_up_after = 5

    def function_and_gradient(self, p):
      objective, gradient = \
        scipy_iterations.scipy_iterations.function_and_gradient(self, p)
      if self.n_iterations >= self.blow_up_after:
        gradient = numpy.full_like(gradient, numpy.nan)
      return objective, gradient

  problem = test_problems.exponential_fit()
  cycles = blowing_up(problem, method='L-BFGS-B', n_max_iterations=500)
  assert cycles.stop_reason is not None and 'finite' in cycles.stop_reason, \
    cycles.stop_reason
  assert cycles.n_iterations >= blowing_up.blow_up_after, cycles.n_iterations
  for value in problem.x:
    assert not math.isnan(value) and not math.isinf(value), list(problem.x)
  # and the objective is still the one belonging to those parameters
  problem.build_up()
  assert not math.isnan(problem.objective()), problem.objective()


def exercise_unsupported_method():
  """ An unknown method shall be reported rather than handed on to scipy. """
  from scitbx.lstbx import scipy_iterations
  problem = test_problems.exponential_fit()
  try:
    scipy_iterations.scipy_iterations(problem, method='no-such-method')
  except RuntimeError as e:
    assert 'no-such-method' in str(e), e
  else:
    raise Exception('an unsupported method went unnoticed')


def run():
  try:
    import scipy.optimize # noqa: F401
  except ImportError:
    print('Skipping tst_scipy_iterations: scipy is not available')
    return
  exercise_gradient_consistency()
  exercise_hessian_product()
  exercise_unsupported_method()
  exercise_non_finite_guard()
  exercise_minimisers()
  exercise_ill_conditioned_problem()
  exercise_preconditioning()
  exercise_interrupt()
  print('OK')


if __name__ == '__main__':
  run()
