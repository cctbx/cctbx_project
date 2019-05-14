from __future__ import absolute_import, division, print_function
from scitbx import minpack
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal

def exercise_interface():
  x = flex.double([1])
  minimizer = minpack.levenberg_marquardt(
    m=1,
    x=x,
    ftol=1.e-16,
    xtol=1.e-16,
    gtol=0,
    factor=1.0e2,
    call_back_after_iteration=True)
  assert minimizer.m == 1
  assert minimizer.ftol == 1.e-16
  assert minimizer.xtol == 1.e-16
  assert minimizer.gtol == 0
  assert minimizer.maxfev == 0
  assert minimizer.factor == 1.0e2
  assert minimizer.requests_fvec()
  assert approx_equal(x, [1])
  minimizer.process_fvec(flex.double([-2.0]))
  assert minimizer.requests_fjac()
  assert approx_equal(x, [1])
  minimizer.process_fjac(flex.double([-1.0]))
  assert minimizer.requests_fvec()
  assert approx_equal(x, [-1])
  minimizer.process_fvec(flex.double([0.0]))
  minimizer.calls_back_after_iteration()
  assert approx_equal(x, [-1])
  minimizer.continue_after_call_back_after_iteration()
  assert minimizer.requests_fjac()
  assert approx_equal(x, [-1])
  minimizer.process_fjac(flex.double([-1.0]))
  assert approx_equal(x, [-1])
  assert minimizer.has_terminated()
  assert minimizer.info == 4
  assert minimizer.nfev == 2
  assert minimizer.njev == 2

def exercise():
  exercise_interface()
  # more comprehensive tests: scitbx/examples/immoptibox_ports.py
  print("OK")

if (__name__ == "__main__"):
  exercise()
