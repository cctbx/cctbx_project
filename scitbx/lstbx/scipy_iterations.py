""" Solve non-linear L.S. problems with the minimisers of scipy.optimize.

The normal equations objects of lstbx already expose everything
scipy.optimize.minimize asks for: build_up() evaluates the objective and its
gradient at the current parameter vector, and normal_matrix_packed_u() is the
Gauss-Newton approximation of the Hessian. This module is the glue.

Unlike the solvers of normal_eqns_solving, which all move along the step
obtained by solving the normal equations, the minimisers driven from here move
the parameters to points of their own choosing. The problem must therefore
implement apply_shifts(), c.f. normal_eqns.non_linear_ls_mixin.

scipy is imported lazily so that importing scitbx.lstbx does not require it.
"""
from __future__ import absolute_import, division, print_function

from scitbx.array_family import flex
from scitbx.lstbx import normal_eqns_solving
import math


class minimisation_stopped(Exception):
  """ Raised to abandon a running minimisation.

  scipy 1.11 lets a callback stop minimize() by raising StopIteration but that
  is too recent to rely on, so this is raised instead and do() catches it.
  Whatever scipy was in the middle of is discarded and the parameters are moved
  back to the last iterate the minimiser accepted.
  """


class method_traits(object):
  """ What scipy.optimize.minimize accepts for one particular method.

  scipy rejects arguments a method cannot consume, and each method spells the
  same tolerance differently, so both have to be spelt out. Supporting a
  further method is a matter of adding one entry to
  scipy_iterations.supported_methods.
  """

  def __init__(self, arguments=(), tolerances=None, max_evaluations=None,
               self_scaling=False):
    """ arguments: which of 'jac', 'hessp' and 'bounds' the method takes.
        tolerances: maps 'f', 'g' and 'x' onto the options the method uses for
          them; those it has no option for are simply left out.
        max_evaluations: the option limiting the number of function
          evaluations, or None if the method has none.
        self_scaling: whether the method works out the size of its own steps,
          either from a trust region it adapts or from a Newton direction. The
          others start from a step of unit norm, which is only a sensible one
          if p is calibrated for it; c.f. calibrate().
    """
    self.arguments = frozenset(arguments)
    self.tolerances = tolerances or {}
    self.max_evaluations = max_evaluations
    self.self_scaling = self_scaling


class scipy_iterations(normal_eqns_solving.iterations):
  """ Iterations driven by scipy.optimize.minimize.

  The minimiser does not work with the parameters x of the problem but with

      p = D (x - x_0)

  where x_0 is the starting point and D the diagonal preconditioner. Working
  with a displacement means the absolute parameter vector is never needed, only
  the increments apply_shifts() takes. The preconditioner matters a great deal:
  Gauss-Newton is invariant under a rescaling of the parameters whereas these
  minimisers are not at all, and the parameters of a crystallographic problem
  differ by orders of magnitude in curvature (fractional coordinates against
  U_iso, say).

  Its magnitude matters as much as its shape, which is what calibrate() is
  about, and a point the model cannot be evaluated at all is reported to the
  minimiser rather than allowed to end the minimisation, which is what
  penalty() is about.

  Note that the histories journaled by the base class are recorded once per
  *evaluation*, of which there are several per iteration.
  """

  method = 'L-BFGS-B'
  preconditioning = 'diagonal' # or None
  calibration = True
  # a run of this many trial points the model cannot be evaluated at is taken
  # to mean the model itself is broken rather than the minimiser overreaching
  max_consecutive_unevaluable = 20
  f_tolerance = None
  g_tolerance = None
  x_tolerance = None
  max_evaluations = None
  scipy_options = None # merged last, so it overrides all of the above
  interrupted = False

  # journaling every trial point costs a vector per evaluation and there are
  # many more evaluations here than there are iterations
  track_step = False
  track_gradient = False

  supported_methods = {
    'CG':           method_traits(('jac',), {'g': 'gtol'}),
    'BFGS':         method_traits(('jac',), {'g': 'gtol'}),
    'L-BFGS-B':     method_traits(('jac', 'bounds'),
                                  {'g': 'gtol', 'f': 'ftol'}, 'maxfun'),
    'TNC':          method_traits(('jac', 'bounds'),
                                  {'g': 'gtol', 'f': 'ftol', 'x': 'xtol'},
                                  'maxfun'),
    'SLSQP':        method_traits(('jac', 'bounds'), {'f': 'ftol'}),
    'Newton-CG':    method_traits(('jac', 'hessp'), {'x': 'xtol'},
                                  self_scaling=True),
    'trust-ncg':    method_traits(('jac', 'hessp'), {'g': 'gtol'},
                                  self_scaling=True),
    'trust-krylov': method_traits(('jac', 'hessp'), {'g': 'gtol'},
                                  self_scaling=True),
    'trust-constr': method_traits(('jac', 'hessp', 'bounds'),
                                  {'g': 'gtol', 'x': 'xtol'},
                                  self_scaling=True),
    'Powell':       method_traits((), {'f': 'ftol', 'x': 'xtol'}, 'maxfev'),
    'Nelder-Mead':  method_traits((), {'f': 'fatol', 'x': 'xatol'}, 'maxfev'),
  }

  def do(self):
    from scipy import optimize
    self.n_iterations = 0
    self.n_evaluations = 0
    self.scipy_result = None
    self.stop_reason = None
    self.objective_scale = 1
    self.unevaluable = None
    self.unevaluable_reason = None
    self.n_unevaluable = 0
    self.consecutive_unevaluable = 0
    traits = self.traits()

    # The first build gives the objective, the gradient and the normal matrix
    # at the starting point; the preconditioner is made out of the latter.
    self.non_linear_ls.build_up()
    self.n_evaluations = 1
    if self.has_gradient_converged_to_zero(): return

    n_parameters = self.non_linear_ls.opposite_of_gradient().size()
    self.parameter_scale = self.preconditioner()
    self.objective_scale = self.calibrate(traits)
    self.x_0 = self.non_linear_ls.parameter_vector()
    self.p_model = flex.double(n_parameters, 0)     # where the parameters are
    self.p_accepted = flex.double(n_parameters, 0)  # last accepted iterate
    self.p_evaluated = flex.double(n_parameters, 0) # where the last build was
    self.p_good = flex.double(n_parameters, 0)      # last one which evaluated
    self.f_good = self.non_linear_ls.objective()*self.objective_scale
    self.evaluated_objective_only = False
    if self.n_max_iterations == 0: return

    keywords = self.minimize_keywords(traits)
    try:
      self.scipy_result = optimize.minimize(
        x0=self.p_model.as_numpy_array(), method=self.method,
        callback=self.on_iteration, **keywords)
    except minimisation_stopped:
      pass
    else:
      self.p_accepted = self.as_flex(self.scipy_result.x)

    # Leave the parameters at the point the minimiser settled on and make the
    # normal equations describe it: the covariance matrix, and with it the
    # standard uncertainties, are worked out from them afterwards.
    self.move_to(self.p_accepted)
    try:
      self.non_linear_ls.build_up()
    except RuntimeError:
      # A minimiser only settles on a point it evaluated, so this is not
      # expected; falling back to one which is known to evaluate is still
      # better than leaving the structure where it cannot be worked with.
      if self.p_accepted.all_eq(self.p_good): raise
      self.stop_reason = "the point the minimiser settled on could not be" \
                         " evaluated; the last one which could was kept"
      self.p_accepted = self.p_good
      self.move_to(self.p_accepted)
      self.non_linear_ls.build_up()

  def __str__(self):
    return "scipy.optimize.minimize, method %s" % self.method

  def traits(self):
    """ The traits of self.method, complaining if it is not a known one. """
    try:
      return self.supported_methods[self.method]
    except KeyError:
      raise RuntimeError(
        "unsupported scipy minimiser '%s'; supported are: %s"
        % (self.method, ', '.join(sorted(self.supported_methods))))

  def preconditioner(self):
    """ The diagonal D such that the minimiser works with p = D (x - x_0).

    D = sqrt(diag(A)) makes the diagonal of the Gauss-Newton Hessian unity.
    """
    if self.preconditioning is None: return None
    if self.preconditioning != 'diagonal':
      raise RuntimeError(
        "unknown preconditioning '%s'" % self.preconditioning)
    diagonal = \
      self.non_linear_ls.normal_matrix_packed_u().matrix_packed_u_diagonal()
    # a parameter the data says nothing about would otherwise scale to infinity
    diagonal.set_selected(diagonal <= 0, 1.)
    return flex.sqrt(diagonal)

  def calibrate(self, traits):
    """ Rescale p so that a step of unit norm in it is the natural step.

    The preconditioner fixes the *shape* of p-space -- the Gauss-Newton Hessian
    has unit diagonal there -- but says nothing about its magnitude, and
    several of these minimisers (L-BFGS-B, TNC, SLSQP) take a step of unit norm
    as their first trial point, having no way to know any better. That step is
    1/||g_p|| times the Gauss-Newton step, so it overshoots worst exactly where
    a refinement is normally started from: an already converged structure,
    whose gradient is small. Overshooting by a factor of ten is enough to drive
    an ADP so far that the structure factors can no longer be worked out at
    all, which is how this surfaces.

    Dividing D by s = ||g_p|| shrinks that first step to the Gauss-Newton one.
    Multiplying the objective by 1/s^2 at the same time is what keeps the units
    honest: the two together leave the unit diagonal alone and make the
    gradient at the starting point unit-norm as well, so the Newton step in p
    is -g_p, of norm 1, which is precisely what those minimisers assume. The
    gradient the g_tolerance is compared against therefore still means what it
    meant.

    A method which scales itself is left alone: it would only have the size of
    its first trial step changed, for no reason, and a trust region converging
    on the right radius from a different one costs iterations.

    Returns the factor to multiply the objective, and its gradient and Hessian
    with it, by; self.parameter_scale is rescaled in place.
    """
    if traits.self_scaling: return 1
    if not self.calibration or self.parameter_scale is None: return 1
    g_p = self.non_linear_ls.opposite_of_gradient()/self.parameter_scale
    s = g_p.norm()
    if s == 0: return 1
    self.parameter_scale = self.parameter_scale/s
    return 1/s**2

  def bounds(self):
    """ Bounds on p for the methods which support them, None for no bounds. """
    return None

  def on_iteration_completion(self):
    """ Called once per iteration accepted by the minimiser; a no-op here. """

  def minimize_keywords(self, traits):
    """ The arguments to pass to scipy.optimize.minimize. """
    options = {'maxiter': self.n_max_iterations}
    for name, value in (('f', self.f_tolerance),
                        ('g', self.g_tolerance),
                        ('x', self.x_tolerance)):
      if value is not None and name in traits.tolerances:
        options[traits.tolerances[name]] = value
    if self.max_evaluations is not None and traits.max_evaluations is not None:
      options[traits.max_evaluations] = self.max_evaluations
    if self.scipy_options is not None:
      options.update(self.scipy_options)
    keywords = {'options': options}
    if 'jac' in traits.arguments:
      keywords['fun'] = self.function_and_gradient
      keywords['jac'] = True
    else:
      keywords['fun'] = self.function
    if 'hessp' in traits.arguments:
      keywords['hessp'] = self.hessian_product
    bounds = self.bounds()
    if bounds is not None:
      if 'bounds' not in traits.arguments:
        raise RuntimeError("'%s' does not support bounds" % self.method)
      keywords['bounds'] = bounds
    return keywords

  def function_and_gradient(self, p):
    self.build_at(p)
    if self.unevaluable is not None:
      objective, gradient = self.penalty(self.as_flex(p))
      return objective, gradient.as_numpy_array()
    gradient = -self.non_linear_ls.opposite_of_gradient()
    if self.parameter_scale is not None:
      gradient /= self.parameter_scale
    gradient *= self.objective_scale
    objective = self.non_linear_ls.objective()*self.objective_scale
    self.check_is_finite(objective, gradient.norm())
    return objective, gradient.as_numpy_array()

  def function(self, p):
    self.build_at(p, objective_only=True)
    if self.unevaluable is not None:
      return self.penalty(self.as_flex(p))[0]
    objective = self.non_linear_ls.objective()*self.objective_scale
    self.check_is_finite(objective)
    return objective

  def penalty(self, p):
    """ What to report at a point where the model cannot be evaluated.

    A minimiser can take the parameters somewhere the model has nothing to say
    at all -- an ADP driven so large that the Debye-Waller factor overflows,
    say -- and there is then no objective to hand back. Letting the failure end
    the minimisation throws away every cycle already done and leaves the
    structure wherever the trial step happened to put it, when all the point
    needs is to be reported as a bad one so that the line search steps back
    from it.

    The surrogate is the quadratic model p-space is preconditioned to look
    like: unit curvature about the last point which did evaluate. It is worse
    there than at that point and its gradient points away from it, which is
    what a line search needs in order to back off.
    """
    d = p - self.p_good
    return self.f_good + 0.5*d.norm()**2, d

  def stop(self, reason):
    """ Abandon the minimisation, keeping the last accepted iterate. """
    self.stop_reason = reason
    raise minimisation_stopped()

  def check_is_finite(self, *values):
    """ Abandon the minimisation rather than let a non-finite value spread.

    The parameters can be taken somewhere the model cannot be evaluated at all
    -- a structure factor overflowing, say -- and every minimiser here would
    then go on to derive its next point from that NaN. It also happens that a
    minimiser divides by a gradient it has itself driven to zero, i.e. on
    success. Stopping gives a worse answer than convergence but a usable one.
    """
    for value in values:
      if math.isnan(value) or math.isinf(value):
        self.stop("the objective or its gradient stopped being finite")

  def hessian_product(self, p, v):
    """ The product of the Gauss-Newton Hessian with v, in units of p.

    The normal matrix is that approximation of the Hessian, so the Newton-type
    minimisers get it for nothing. It includes the floating origin restraint,
    which contributes nothing to the gradient but makes the directions along
    which the objective is flat well determined.
    """
    self.build_at(p)
    u = self.as_flex(v)
    if self.unevaluable is not None:
      # the unit curvature the penalty above is built on
      return u.as_numpy_array()
    if self.parameter_scale is not None:
      u /= self.parameter_scale
    u.reshape(flex.grid(1, u.size()))
    product = u.matrix_multiply_packed_u(
      self.non_linear_ls.normal_matrix_packed_u())
    product.reshape(flex.grid(product.size()))
    if self.parameter_scale is not None:
      product /= self.parameter_scale
    product *= self.objective_scale
    return product.as_numpy_array()

  def on_iteration(self, p, *args):
    """ Called by scipy once per accepted iteration.

    trust-constr passes the optimisation state as a second argument, the other
    methods pass the iterate alone.
    """
    self.n_iterations += 1
    self.p_accepted = self.as_flex(p)
    self.on_iteration_completion()
    if self.interrupted: self.stop("the minimisation was interrupted")

  def build_at(self, p, objective_only=False):
    """ Make the normal equations describe the point p, rebuilding if needed.

    scipy asks for the objective, the gradient and Hessian products at the same
    point several times over, and each build is a whole pass over the data.
    """
    p = self.as_flex(p)
    if (self.p_evaluated is not None and p.all_eq(self.p_evaluated)
        and (objective_only or not self.evaluated_objective_only)):
      return
    if math.isnan(p.norm()) or math.isinf(p.norm()):
      # Moving there would leave the parameters non-finite and, the moves being
      # increments of one another, past recovery: the last accepted iterate
      # could no longer be reached by shifting back from a NaN.
      self.stop("the minimiser proposed a point which is not finite")
    if (self.max_evaluations is not None
        and self.n_evaluations >= self.max_evaluations):
      self.stop("the limit of %i evaluations was reached"
                % self.max_evaluations)
    self.move_to(p)
    try:
      self.non_linear_ls.build_up(objective_only=objective_only)
    except RuntimeError as e:
      self.consecutive_unevaluable += 1
      if self.consecutive_unevaluable > self.max_consecutive_unevaluable:
        # every trial point failing is the model being broken, not a line
        # search overreaching, and there is nothing useful to report back
        raise
      self.n_unevaluable += 1
      self.unevaluable = e
      self.unevaluable_reason = str(e)  # kept once self.unevaluable is cleared
    else:
      self.consecutive_unevaluable = 0
      self.unevaluable = None
      self.p_good = p
      self.f_good = self.non_linear_ls.objective()*self.objective_scale
    self.n_evaluations += 1
    self.p_evaluated = p
    self.evaluated_objective_only = objective_only

  def invalidate_evaluation(self):
    """ Have the next evaluation rebuild the normal equations.

    Whatever consumes them has to say so: solving them for a covariance matrix,
    for instance, takes the normal matrix away, and the cache would otherwise
    go on believing it is there to be handed out.
    """
    self.p_evaluated = None

  def move_to(self, p):
    """ Move the parameters to where the minimiser wants them. """
    shifts = p - self.p_model
    if shifts.all_eq(0): return
    if self.parameter_scale is not None:
      shifts /= self.parameter_scale
    self.non_linear_ls.apply_shifts(shifts)
    self.p_model = p.deep_copy()
    x = self.non_linear_ls.parameter_vector()
    if x is not None and self.x_0 is not None:
      # Parameters get clamped as the shifts are applied -- a negative U_iso,
      # say -- so where they end up is not necessarily where the minimiser
      # asked them to go. Taking the request at face value would leave every
      # later shift off by the difference, in particular the ones undoing a
      # trial step the minimiser went on to reject.
      self.p_model = x - self.x_0
      if self.parameter_scale is not None:
        self.p_model *= self.parameter_scale

  @staticmethod
  def as_flex(a):
    """ scipy hands out numpy arrays, which flex takes directly. """
    if isinstance(a, flex.double): return a.deep_copy()
    return flex.double(a)
