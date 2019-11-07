from __future__ import absolute_import, division, print_function
import scitbx.math
from scitbx.array_family import flex
from libtbx.test_utils import Exception_expected, approx_equal
from libtbx.utils import format_cpu_times
import time
import sys
from six.moves import range

def exercise_with_random_arguments(n_arguments, n_iterations):
  mt = flex.mersenne_twister(seed=0)
  d = mt.random_double(size=n_arguments)*100-50
  f = d.as_float()
  print("showing wall clock times!")
  t0 = time.time()
  for i_iteration in range(n_iterations):
    jef = scitbx.math.jacks_expf(array_of_float=f)
  print("jacks_expf(): %.2f s" % (time.time()-t0))
  for i_iteration in range(n_iterations):
    sef = flex.exp(f)
  print("std::exp(float): %.2f s" % (time.time()-t0))
  t0 = time.time()
  for i_iteration in range(n_iterations):
    sed = flex.exp(d)
  print("std::exp(double): %.2f s" % (time.time()-t0))
  assert sed.all_gt(0)
  for ef in [jef, sef]:
    max_rel_err = flex.max(flex.abs((ef.as_double() - sed) / sed))
    assert approx_equal(max_rel_err, 0, eps=1e-5)

def run(args):
  assert args in [[], ["--comprehensive"]]
  if (len(args) == 0):
    comprehensive = False
  else:
    comprehensive = True
  if (not comprehensive):
    n_arguments = 100000
    n_iterations = 10
    exponents = [-127, -1, 0, 1, 5, 6]
    mantissa_step_size = 100000
    j_sample = 300000
  else:
    n_arguments = 300000
    n_iterations = 20
    exponents = range(-127, 7)
    mantissa_step_size = 1
    j_sample = 300000
  exercise_with_random_arguments(
    n_arguments=n_arguments,
    n_iterations=n_iterations)
  for negative_sign in [False, True]:
    for exponent in exponents:
      if (comprehensive): print(exponent)
      try:
        r = scitbx.math.exercise_jacks_expf(
          negative_sign=negative_sign,
          exponent=exponent,
          mantissa_step_size=mantissa_step_size,
          j_sample=j_sample)
        if (comprehensive): print(list(r))
      except RuntimeError as e:
        assert str(e) == "jacks_expf(): function argument out of range."
        assert not negative_sign
        break
    else:
      if (not negative_sign):
        raise Exception_expected
  #
  print(format_cpu_times())

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
