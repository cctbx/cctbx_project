from __future__ import absolute_import, division, print_function
import math
from scitbx.array_family import flex
from scitbx import lbfgs
from six.moves import range

def exercise_drop_convergence_test():
  c = lbfgs.drop_convergence_test()
  assert c.n_test_points() > 0
  assert c.max_drop_eps() > 0
  assert c.iteration_coefficient() > 0
  assert c.objective_function_values().size() == 0
  assert c.max_drop() == 0
  c = lbfgs.drop_convergence_test(n_test_points=6)
  c = lbfgs.drop_convergence_test(n_test_points=6, max_drop_eps=1.e-3)
  c = lbfgs.drop_convergence_test(
    n_test_points=6, max_drop_eps=1.e-3, iteration_coefficient=3)
  assert c.n_test_points() == 6
  assert c.max_drop_eps() == 1.e-3
  assert c.iteration_coefficient() == 3
  assert c.objective_function_values().size() == 0
  assert c.max_drop() == 0
  for n_test_points in (2, 7):
    c = lbfgs.drop_convergence_test(n_test_points, 1.e-3)
    assert c.n_test_points() == n_test_points
    converged = []
    for x in range(10):
      converged.append(c(math.exp(-x)))
    c.objective_function_values().size() == 10
    if (n_test_points == 2):
      assert converged == [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
    else:
      assert converged == [0, 0, 0, 0, 0, 0, 0, 0, 1, 1]

def exercise_minimization(verbose):
  n = 100
  x = flex.double(n)
  g = flex.double(n)
  for j in range(0, n, 2):
    x[j] = -1.2
    x[j+1] = 1.
  minimizer = lbfgs.minimizer(n)
  is_converged = lbfgs.traditional_convergence_test(n)
  if (verbose):
    print("n: ", minimizer.n())
    print("m: ", minimizer.m())
    print("xtol: ", minimizer.xtol())
  while 1:
    f = 0.
    for j in range(0, n, 2):
      t1 = 1.e0 - x[j]
      t2 = 1.e1 * (x[j+1] - x[j] * x[j])
      g[j+1] = 2.e1 * t2
      g[j] = -2.e0 * (x[j] * g[j+1] + t1)
      f = f + t1 * t1 + t2 * t2
    if (minimizer.run(x, f, g)): continue
    if (verbose):
      print("f:", f, "gnorm:", minimizer.euclidean_norm(g))
      print(minimizer.iter(), minimizer.nfun(), minimizer.stp())
    if (is_converged(x, g)): break
    if (minimizer.nfun() > 2000): break
    assert minimizer.run(x, f, g)
  assert f < 1.e-12
  assert minimizer.euclidean_norm(g) < 1.e-4
  assert minimizer.iter() < 40
  assert minimizer.nfun() < 50
  assert abs(minimizer.stp() - 1) < 1.e-6

def exercise_first_step_with_diagonal():
  # With the inverse diagonal of the Hessian supplied, the search
  # direction of the first iteration is the Newton step, so the first
  # trial point of a quadratic is its minimum. Without a diagonal the
  # first trial is a step of length 1 along -g (Nocedal's heuristic).
  c = flex.double([1., 10., 100.])
  x = flex.double([1., 2., 3.])
  g = c * x
  f = 0.5 * flex.sum(c * x * x)
  minimizer = lbfgs.minimizer(3)
  assert minimizer.run(x, f, g, 1. / c)
  assert x.all_approx_equal(0., 1.e-12), list(x)
  x = flex.double([1., 2., 3.])
  g = c * x
  minimizer = lbfgs.minimizer(3)
  assert minimizer.run(x, f, g)
  expected = flex.double([1., 2., 3.]) - g / minimizer.euclidean_norm(g)
  assert x.all_approx_equal(expected, 1.e-12), list(x)

def run():
  import os, sys
  Endless = "--Endless" in sys.argv
  verbose = "--Verbose" in sys.argv and not Endless
  exercise_drop_convergence_test()
  exercise_first_step_with_diagonal()
  while 1:
    exercise_minimization(verbose)
    if (not Endless): break
  t = os.times()
  print("OK")

if (__name__ == "__main__"):
  run()
