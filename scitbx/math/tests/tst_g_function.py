from __future__ import absolute_import, division, print_function
from scitbx import math as m
from libtbx.test_utils import approx_equal
import math

def G(x):
  if(x != 0):
    return 3*(math.sin(x)-x*math.cos(x))/x**3
  else:
    return 1

def exercise():
  x = 0
  while x < 4.45: #undefined for larger arguments
    gf = m.GfuncOfRSsqr_approx((x/(2.*math.pi))**2)
    assert approx_equal(gf, G(x))
    x += 0.01

if (__name__ == "__main__"):
  exercise()
