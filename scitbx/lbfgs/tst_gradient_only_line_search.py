"""
scitbx.lbfgs gradient-only mode: the line search accepts a step only on the
strong curvature condition, brackets an overshoot, expands an undershoot,
never evaluates the same point twice, takes a first trial step of
first_step_max_shift, and stores only pairs with positive curvature. The
targets return f for the test's own checks; the minimizer sees gradients
only.
"""
import math
import scitbx.lbfgs
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from scitbx.lbfgs.tst_max_shift import largest_shift

H = flex.double([1e4, 1, 1, 1e4, 1, 1])
C = flex.double([1, 2, 3, 4, 5, 6])

def ill_quadratic(x):
  d = x - C
  return flex.sum(H * d * d), 2 * H * d

def rosenbrock3(x):
  x0, x1, x2 = x
  f = 100 * (x1 - x0**2)**2 + (1 - x0)**2 + (x2 - 1)**2
  g = flex.double([
    -400 * x0 * (x1 - x0**2) - 2 * (1 - x0),
    200 * (x1 - x0**2),
    2 * (x2 - 1)])
  return f, g

def gaussian_wells(x):
  s2 = 2 * 0.3**2
  f = 0
  g = flex.double(3, 0)
  for c in [(0, 0, 0), (1.2, 0, 0)]:
    d = x - flex.double(c)
    e = math.exp(-flex.sum(d * d) / s2)
    f -= e
    g += 2 * d * e / s2
  return f, g

def abs_target(x):
  g = flex.double([1 if v > 0 else -1 for v in x])
  return flex.sum(flex.abs(x)), g

def two_wells(x):
  f = 0
  g = flex.double(x.size(), 0)
  for i, v in enumerate(x):
    if(v < 1):
      f += v * v
      g[i] = 2 * v
    else:
      f += (v - 2)**2
      g[i] = 2 * (v - 2)
  return f, g

class evaluator(object):
  def __init__(self, target, x0):
    self.target = target
    self.x = flex.double(x0)
    self.x_history = []
    self.f_history = []
    self.iterates = []
    self.g = None

  def compute_functional_and_gradients(self):
    f, self.g = self.target(self.x)
    self.x_history.append(self.x.deep_copy())
    self.f_history.append(f)
    if(len(self.iterates) == 0):
      self.iterates.append((self.x.deep_copy(), self.g.deep_copy()))
    return f, self.g

  def callback_after_step(self, minimizer):
    self.iterates.append((self.x.deep_copy(), self.g.deep_copy()))
    return False

def run_go(target, x0, line_search=True, max_iterations=500, **core):
  e = evaluator(target, x0)
  m = scitbx.lbfgs.run(target_evaluator=e, gradient_only=True,
    line_search=line_search,
    core_params=scitbx.lbfgs.core_parameters(**core),
    termination_params=scitbx.lbfgs.termination_parameters(
      max_iterations=max_iterations))
  return e, m

def exercise_converges():
  for target, x0, xmin in [(ill_quadratic, [0] * 6, list(C)),
                           (rosenbrock3, [-1.2, 1, 0], [1, 1, 1]),
                           (gaussian_wells, [0.5, 0.1, 0], [0, 0, 0])]:
    e, m = run_go(target, x0, max_shift=0.2, first_step_max_shift=0.05)
    assert approx_equal(e.x, xmin, 1e-3), (target.__name__, list(e.x))
    assert e.f_history[-1] < e.f_history[0]

def exercise_no_duplicate_evaluation():
  for target, x0 in [(ill_quadratic, [0] * 6), (gaussian_wells, [0.5, 0.1, 0])]:
    e, m = run_go(target, x0, max_shift=0.2, first_step_max_shift=0.05)
    for a, b in zip(e.x_history, e.x_history[1:]):
      assert (a - b).norm() > 0, "same point evaluated twice"

def exercise_curvature_pairs():
  e, m = run_go(ill_quadratic, [0] * 6, max_shift=0.2,
    first_step_max_shift=0.05)
  assert len(e.iterates) > 2
  for (x0, g0), (x1, g1) in zip(e.iterates, e.iterates[1:]):
    s = x1 - x0
    y = g1 - g0
    assert flex.sum(s * y) > 0

def exercise_first_step():
  e, m = run_go(ill_quadratic, [0] * 6, max_shift=0.2,
    first_step_max_shift=0.05)
  first = e.x_history[1] - e.x_history[0]
  assert approx_equal(largest_shift(first), 0.05, 1e-9), largest_shift(first)

def exercise_exhaustion_and_guard():
  e, m = run_go(abs_target, [0.3, -0.2, 0.1], max_shift=0.2,
    first_step_max_shift=0.05)
  assert m.n_pairs() <= 5
  assert e.f_history[-1] <= e.f_history[0]

def exercise_non_descent_reset():
  e, m = run_go(two_wells, [0.9, 0.9, 0.9], max_shift=0.5,
    first_step_max_shift=0.05)
  assert e.f_history[-1] < e.f_history[0]

def exercise_zero_gradient_start():
  e, m = run_go(ill_quadratic, list(C), max_shift=0.2,
    first_step_max_shift=0.05)
  assert len(e.x_history) == 1, len(e.x_history)

def exercise_no_line_search_flag():
  e, m = run_go(ill_quadratic, [0] * 6, line_search=False, max_shift=0.2,
    first_step_max_shift=0.05)
  # one evaluation per iteration, plus the start
  assert len(e.x_history) == len(e.iterates), (
    len(e.x_history), len(e.iterates))

def flat_basin(x):
  """Zero gradient inside |x_i| < 1, quadratic walls outside."""
  f = 0
  g = flex.double(x.size(), 0)
  for i, v in enumerate(x):
    if(v > 1):
      f += (v - 1)**2
      g[i] = 2 * (v - 1)
    elif(v < -1):
      f += (v + 1)**2
      g[i] = 2 * (v + 1)
  return f, g

def exercise_exhaustion_progress():
  """When no trial point reduces |d| (kink of a non-smooth target) the
  search keeps a non-zero step instead of freezing x: the minimum is
  reached."""
  e, m = run_go(abs_target, [0.3, 0.3, 0.3], max_iterations=20,
    max_shift=0.2, first_step_max_shift=0.05)
  assert flex.max(flex.abs(e.x)) < 1e-6, list(e.x)

def exercise_zero_gradient_after_step():
  """A zero gradient after an accepted step converges, no internal error,
  even when min_iterations keeps the Python loop from stopping first."""
  e = evaluator(flat_basin, [1.03, 1.03, 1.03])
  scitbx.lbfgs.run(target_evaluator=e, gradient_only=True,
    core_params=scitbx.lbfgs.core_parameters(max_shift=0.2,
      first_step_max_shift=0.05),
    termination_params=scitbx.lbfgs.termination_parameters(
      min_iterations=3, max_iterations=20))
  assert flex.max(flex.abs(e.x)) <= 1 + 1e-12, list(e.x)

def concave(x):
  """-sum x^2: every curvature pair has s.y < 0 and is skipped, so the
  second iteration takes the steepest-descent reset path."""
  return -flex.sum(x * x), -2 * x

class diag_evaluator(evaluator):
  """Supplies a fixed diagonal (inverse Hessian estimate) every iteration."""
  diag_mode = "always"
  def __init__(self, target, x0, diag):
    evaluator.__init__(self, target, x0)
    self.diag = flex.double(diag)
  def compute_functional_gradients_diag(self):
    f, g = self.compute_functional_and_gradients()
    return f, g, self.diag

def exercise_diag_after_reset():
  """A user diagonal shapes the steepest-descent direction after a reset
  as it does in the first iteration: with diag (1, 1e-6, 1e-6) every
  accepted step moves x1 and x2 by at most 1e-3 of x0's move."""
  e = diag_evaluator(concave, [1, 1, 1], diag=[1, 1e-6, 1e-6])
  scitbx.lbfgs.run(target_evaluator=e, gradient_only=True,
    core_params=scitbx.lbfgs.core_parameters(max_shift=0.2,
      first_step_max_shift=0.05),
    termination_params=scitbx.lbfgs.termination_parameters(
      max_iterations=3))
  assert len(e.iterates) == 4, len(e.iterates)
  for (x0, g0), (x1, g1) in zip(e.iterates, e.iterates[1:]):
    s = x1 - x0
    assert abs(s[0]) > 0
    assert max(abs(s[1]), abs(s[2])) <= 1e-3 * abs(s[0]), list(s)

def exercise_function_mode_unchanged():
  e1 = evaluator(rosenbrock3, [-1.2, 1.0, 0.5])
  scitbx.lbfgs.run(target_evaluator=e1)
  e2 = evaluator(rosenbrock3, [-1.2, 1.0, 0.5])
  scitbx.lbfgs.run(target_evaluator=e2,
    core_params=scitbx.lbfgs.core_parameters(first_step_max_shift=0.05))
  assert [list(x) for x in e1.x_history] == [list(x) for x in e2.x_history]

def run():
  exercise_converges()
  exercise_no_duplicate_evaluation()
  exercise_curvature_pairs()
  exercise_first_step()
  exercise_exhaustion_and_guard()
  exercise_non_descent_reset()
  exercise_zero_gradient_start()
  exercise_no_line_search_flag()
  exercise_exhaustion_progress()
  exercise_zero_gradient_after_step()
  exercise_diag_after_reset()
  exercise_function_mode_unchanged()
  print("OK")

if (__name__ == "__main__"):
  run()
