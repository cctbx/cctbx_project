"""
scitbx.lbfgs core_parameters.max_shift: the largest shift of any atom
(consecutive coordinate triplet; single coordinates when n is not a multiple
of 3) in one line search is bounded, in function mode as in gradient-only
mode, and max_shift=0 (the default) leaves the iterates unchanged.
"""
import scitbx.lbfgs
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal

def largest_shift(d):
  n = d.size()
  if(n % 3 == 0):
    return flex.max(flex.sqrt(flex.vec3_double(d).dot()))
  return flex.max(flex.abs(d))

def quadratic(x):
  c = flex.double([10, 0, 0, 0, 0, 0])
  return flex.sum((x - c) * (x - c)), 2 * (x - c)

def quadratic4(x):
  c = flex.double([3, -2, 0, 1])
  return flex.sum((x - c) * (x - c)), 2 * (x - c)

def rosenbrock3(x):
  x0, x1, x2 = x
  f = 100 * (x1 - x0**2)**2 + (1 - x0)**2 + (x2 - 1)**2
  g = flex.double([
    -400 * x0 * (x1 - x0**2) - 2 * (1 - x0),
    200 * (x1 - x0**2),
    2 * (x2 - 1)])
  return f, g

class evaluator(object):
  """Records every evaluation (x, f) and the x of every accepted iteration."""
  def __init__(self, target, x0):
    self.target = target
    self.x = flex.double(x0)
    self.history = []
    self.iterates = [self.x.deep_copy()]

  def compute_functional_and_gradients(self):
    f, g = self.target(self.x)
    self.history.append((list(self.x), f))
    return f, g

  def callback_after_step(self, minimizer):
    self.iterates.append(self.x.deep_copy())
    return False

  def steps(self):
    return [largest_shift(b - a)
      for a, b in zip(self.iterates, self.iterates[1:])]

def exercise_cap_function_mode():
  e = evaluator(target=quadratic, x0=[0] * 6)
  scitbx.lbfgs.run(target_evaluator=e,
    core_params=scitbx.lbfgs.core_parameters(max_shift=0.5))
  assert len(e.steps()) > 1
  assert max(e.steps()) <= 0.5 + 1e-12, max(e.steps())
  assert approx_equal(e.x, [10, 0, 0, 0, 0, 0], 1e-4)

def exercise_cap_per_coordinate():
  e = evaluator(target=quadratic4, x0=[0] * 4)
  scitbx.lbfgs.run(target_evaluator=e,
    core_params=scitbx.lbfgs.core_parameters(max_shift=0.1))
  assert len(e.steps()) > 1
  assert max(e.steps()) <= 0.1 + 1e-12, max(e.steps())
  assert approx_equal(e.x, [3, -2, 0, 1], 1e-4)

def exercise_no_cap_identical():
  e1 = evaluator(target=rosenbrock3, x0=[-1.2, 1.0, 0.5])
  scitbx.lbfgs.run(target_evaluator=e1)
  e2 = evaluator(target=rosenbrock3, x0=[-1.2, 1.0, 0.5])
  scitbx.lbfgs.run(target_evaluator=e2,
    core_params=scitbx.lbfgs.core_parameters(max_shift=0))
  assert e1.history == e2.history
  assert approx_equal(e1.x, [1, 1, 1], 1e-4)

def exercise_accessors():
  mm = scitbx.lbfgs.ext.minimizer(6, 5, 20, 0.9, 1e-16, 1e-20, 1e20, 0.5, 0.05)
  assert mm.max_shift() == 0.5
  assert mm.first_step_max_shift() == 0.05
  mm = scitbx.lbfgs.ext.minimizer(6)
  assert mm.max_shift() == 0
  assert mm.first_step_max_shift() == 0

def run():
  exercise_cap_function_mode()
  exercise_cap_per_coordinate()
  exercise_no_cap_identical()
  exercise_accessors()
  print("OK")

if (__name__ == "__main__"):
  run()
