from __future__ import division
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
  while x < 5:
    gf = math.GfuncOfRSsqr(x)
    print "%6.4f %6.4f %6.4f"%(x, gf, G(x))
    x += 0.01

if (__name__ == "__main__"):
  exercise()
