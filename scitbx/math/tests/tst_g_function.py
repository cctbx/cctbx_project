from __future__ import division
from scitbx import math
from stdlib import math as smath

def G(x):
  if(x != 0):
    return 3*(smath.sin(x)-x*smath.cos(x))/x**3
  else:
    return 1

def exercise():
  x = 0
  while x < 20:
    gf = math.GfuncOfRSsqr_approx(x)
    print "%6.4f %6.4f %6.4f"%(x, gf, G(x))
    x += 0.01


if (__name__ == "__main__"):
  exercise()
