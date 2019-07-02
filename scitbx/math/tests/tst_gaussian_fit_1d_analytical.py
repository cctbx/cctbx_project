from __future__ import absolute_import, division, print_function
import scitbx.math
from scitbx.array_family import flex
import math
from libtbx.test_utils import approx_equal
from six.moves import range

def exercise_1():
  x = flex.double()
  y = flex.double()
  for x_ in range(0, 10):
    x_ = x_/10.
    x.append(x_)
    y.append( 2.*math.exp(-3.*x_**2) )
  r = scitbx.math.gaussian_fit_1d_analytical(x = x, y = y)
  assert approx_equal(r.a, 2., 1.e-6)
  assert approx_equal(r.b, 3., 1.e-6)

def exercise_2():
  x = flex.double()
  z = flex.double()
  for i in range(1,11):
    x.append(i)
  for i in range(11,21):
    z.append(i)
  y = z * 2 * flex.exp(-5.*x*x)
  r = scitbx.math.gaussian_fit_1d_analytical(x = x, y = y, z = z)
  assert approx_equal(r.a, 2., 1.e-6)
  assert approx_equal(r.b, 5., 1.e-6)

if (__name__ == "__main__"):
  exercise_1()
  exercise_2()
  print("OK")
