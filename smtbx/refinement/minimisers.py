""" Refine a crystal structure with the minimisers of scipy.optimize.

scitbx.lstbx.scipy_iterations does the work; what is added here is the handful
of things that only make sense for a crystallographic L.S. problem: tolerances
of a size the normalised crystallographic objective actually reaches, a report
of how the refinement went, and a warning when a minimiser which needs room to
work is given the four cycles habitually asked of Gauss-Newton.

The objective and its gradient come straight out of the normal equations, so
there is nothing crystallographic about the minimisation itself:

  - crystallographic_ls.objective() is L(K*(x), x), the residual with the
    overall scale factor already optimised away, normalised by sum w Fo^4 and
    including the restraints;
  - its gradient is what the normal equations call the right hand side, which
    accounts for the scale factor having been optimised away;
  - the normal matrix is the Gauss-Newton approximation of the Hessian, which
    the Newton-type minimisers are handed as Hessian-vector products.

A caveat worth knowing about, since it is the one thing here that is not exact.
The SHELX weighting schemes make the weights a function of Fc and of the
overall scale factor, and the gradient above treats them as constants; on top
of that, build_up() hands the weighting scheme the scale factor optimised at
the point evaluated before it. So with those schemes the objective is not quite
a function of the parameters alone, and its gradient is not quite its gradient.
Gauss-Newton never notices, walking towards the minimum and never back, whereas
a line search probing a point and coming back does.

tst_scipy_minimisers.py measures both, and on a small structure with noisy data
finds, for the usual a = 0.1:

  - the gradient out by some 5% of its norm against finite differences;
  - the objective at one and the same point out by some 1% after a detour of
    0.05 in four parameters, recovering to 6e-4 after one further build.

With unit, sigma or stl weights, or SHELX weights with a = b = 0, the weights
depend on neither Fc nor the scale factor and both effects vanish. The gradient
then agrees with finite differences to 2e-10, and every method reaches exactly
the objective Gauss-Newton reaches -- which is the useful thing to know: none
of what follows is a shortcoming of the minimisers.

Where a is not zero they all stop a little short of the Gauss-Newton objective
instead, by an amount which differs from method to method.

Newton-CG is the default because it is handed the Gauss-Newton normal matrix as
Hessian-vector products, which is exactly the information Gauss-Newton itself
uses, and on a crystal structure it reaches the same minimum in a comparable
number of cycles. The scitbx test problems rank the methods differently, but
they are far worse conditioned than a crystallographic refinement and are not a
good guide to this one.

CG and L-BFGS-B use gradients alone and converge more slowly. They are worth
having for problems where the normal matrix is unhelpful, but they are not
substitutes for Gauss-Newton.
"""
from __future__ import absolute_import, division, print_function

from scitbx.lstbx import scipy_iterations
import sys


class crystallographic_scipy_iterations(scipy_iterations.scipy_iterations):
  """ scipy_iterations with defaults suited to a crystal structure refinement.

  The tolerances are relative to an objective which is normalised to lie
  between 0 and 1 and which sits around wR2^2/2, i.e. some 1e-3 for a refined
  structure; asking a minimiser for much more than g_tolerance=1e-8 on such a
  quantity is asking for more than the data support.
  """

  method = 'Newton-CG'
  g_tolerance = 1e-8
  f_tolerance = 1e-12
  # A gradient method makes far less progress per iteration than Gauss-Newton
  # does, so the habitual four cycles are nowhere near enough.
  n_max_iterations = 50
  verbose = False
  log = None

  # Below this many iterations, a method which is not self-scaling is being
  # asked to converge on a budget it cannot converge on, and says so.
  iterations_worth_warning_about = 10
  self_scaling_methods = ('trust-ncg', 'trust-krylov', 'trust-constr',
                          'Newton-CG')

  def do(self):
    if (self.n_max_iterations
        and self.n_max_iterations < self.iterations_worth_warning_about
        and self.method not in self.self_scaling_methods):
      print("Warning: %s is being given only %i iterations; it typically needs"
            " more than Gauss-Newton, not fewer"
            % (self.method, self.n_max_iterations), file=self.out())
    scipy_iterations.scipy_iterations.do(self)
    if self.verbose: self.show_summary()

  def out(self):
    return self.log if self.log is not None else sys.stdout

  def on_iteration_completion(self):
    """ Report on the point last evaluated.

    Which is the iterate just accepted for every method whose line search ends
    on the point it accepts -- all of those used in practice -- and in any case
    the point all three figures below describe, since they come from the same
    build.
    """
    if not self.verbose: return
    normal_eqns = self.non_linear_ls.actual
    print("  %4i  %12.8f  %8.4f  %8.4f  %6d"
          % (self.n_iterations, normal_eqns.objective(), normal_eqns.wR2(),
             normal_eqns.r1_factor()[0], self.n_evaluations), file=self.out())

  def show_summary(self):
    log = self.out()
    print("%s: %i iterations, %i evaluations of the objective"
          % (self, self.n_iterations, self.n_evaluations), file=log)
    if self.stop_reason is not None:
      print("  stopped early because %s" % self.stop_reason, file=log)
    elif self.scipy_result is not None:
      print("  %s" % self.scipy_result.message, file=log)


def shift_over_su(normal_eqns, shifts):
  """ The shifts of the independent parameters in units of their s.u.

  Gauss-Newton gets this from the step it is about to take; a scipy minimiser
  has no such step, so it is worked out from the shifts actually applied. The
  normal equations have to have been built at the point the shifts start from,
  and solving them here is what makes the covariance matrix available.
  """
  from scitbx.array_family import flex
  variances = normal_eqns.covariance_matrix().matrix_packed_u_diagonal()
  variances.set_selected(variances <= 0, 1.)
  return shifts/flex.sqrt(variances)
