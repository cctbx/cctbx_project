from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from scitbx import lbfgsb
import scitbx.minimizers
import math

"""
Termination controls of scitbx.lbfgsb.run and scitbx.minimizers.lbfgs
(mode lbfgsb): factr, pgtol, and what max_iterations counts.
"""

class rosenbrock(object):
  """
  Generalized Rosenbrock function; optionally with noise of fixed size in the
  function value (as a target computed with limited accuracy).
  """
  def __init__(self, n=20, noise=0):
    self.n = n
    self.noise = noise
    self.x = flex.double(n, -1.2)
    self.n_calls = 0
    self.bound_flags = flex.int(n, 2)
    self.lower_bound = flex.double(n, -5)
    self.upper_bound = flex.double(n,  5)

  def _tg(self):
    x = self.x
    t = 0
    g = flex.double(self.n, 0)
    for i in range(self.n-1):
      a = x[i+1] - x[i]**2
      b = 1 - x[i]
      t += 100*a*a + b*b
      g[i]   += -400*a*x[i] - 2*b
      g[i+1] += 200*a
    if(self.noise):
      t += self.noise*math.sin(1.e+6*t)
    return t, g

  # interface of scitbx.lbfgsb.run
  def compute_functional_and_gradients(self):
    self.n_calls += 1
    t, g = self._tg()
    return self.x, t, g

  # interface of scitbx.minimizers.lbfgs
  def update(self, x):
    self.x = x
    self.n_calls += 1
    self.t, self.g = self._tg()

  def target(self): return self.t

  def gradients(self): return self.g

def run_low_level(**kw):
  o = rosenbrock(noise = kw.pop("noise", 0))
  m = lbfgsb.run(
    target_evaluator = o,
    bound_flags      = o.bound_flags,
    lower_bound      = o.lower_bound,
    upper_bound      = o.upper_bound,
    n                = o.n,
    **kw)
  return o, m

def exercise_max_iterations():
  """
  max_iterations limits function evaluations by default, iterations with
  count_iterations=True.
  """
  o, m = run_low_level()
  assert flex.max(flex.abs(o.x - 1)) < 1.e-3 # converged without a limit
  n_iter_full, n_calls_full = m.n_iteration(), o.n_calls
  assert n_calls_full > n_iter_full > 30
  limit = 20
  o, m = run_low_level(max_iterations = limit)
  # stops at the end of the iteration in which the limit on calls is passed
  assert limit < o.n_calls <= limit + 22, o.n_calls
  assert m.n_iteration() < limit
  o, m = run_low_level(max_iterations = limit, count_iterations = True)
  assert m.n_iteration() == limit, m.n_iteration()
  assert o.n_calls > limit

def exercise_factr_pgtol():
  """
  Larger factr or pgtol: earlier stop. With noise in the function value the
  default factr keeps the minimizer going long after it stopped making
  progress; a factr matched to the noise stops it.
  """
  o_default, m = run_low_level()
  o_factr, m = run_low_level(factr = 1.e+12)
  o_pgtol, m = run_low_level(pgtol = 1.)
  assert o_factr.n_calls < o_default.n_calls
  assert o_pgtol.n_calls < o_default.n_calls
  assert flex.max(flex.abs(o_default.x - 1)) < 1.e-3
  # noise 1e-6 in a target of order 1..1e3
  o1, m = run_low_level(noise = 1.e-6, max_iterations = 2000)
  o2, m = run_low_level(noise = 1.e-6, max_iterations = 2000, factr = 1.e+10)
  assert o2.n_calls < o1.n_calls, [o2.n_calls, o1.n_calls]
  assert flex.max(flex.abs(o2.x - 1)) < 0.05

def exercise_minimizers_interface():
  """
  scitbx.minimizers.lbfgs passes factr and pgtol on (lbfgsb mode).
  """
  n_calls = []
  for kw in [dict(), dict(factr = 1.e+12), dict(pgtol = 1.)]:
    o = rosenbrock()
    scitbx.minimizers.lbfgs(mode = "lbfgsb", calculator = o, **kw)
    n_calls.append(o.n_calls)
  assert n_calls[1] < n_calls[0] and n_calls[2] < n_calls[0], n_calls

if(__name__ == "__main__"):
  exercise_max_iterations()
  exercise_factr_pgtol()
  exercise_minimizers_interface()
  print("OK")
