import scitbx.math
from scitbx.array_family import flex
import math
from libtbx.test_utils import approx_equal

def exercise():
  x = flex.double()
  y = flex.double()
  for x_ in xrange(0, 10):
    x_ = x_/10.
    x.append(x_)
    y.append( 2.*math.exp(-3.*x_**2) )
  r = scitbx.math.gaussian_fit_1d_analytical(x = x, y = y)
  assert approx_equal(r.a, 2., 1.e-6)
  assert approx_equal(r.b, 3., 1.e-6)

if (__name__ == "__main__"):
  exercise()
  print "OK"
