import Numeric
from fortran_lbfgs import lbfgs

def run():
  n = 100
  m = 5
  x = Numeric.array(Numeric.arange(n), Numeric.Float64)
  g = Numeric.array(Numeric.arange(n), Numeric.Float64)
  diag = Numeric.array(Numeric.arange(n), Numeric.Float64)
  iprint = [1, 0]
  diagco = 0
  eps = 1.0e-5
  xtol = 1.0e-16
  for j in xrange(0, n, 2):
    x[j] = -1.2
    x[j+1] = 1.
  size_w = n*(2*m+1)+2*m
  w = Numeric.array(Numeric.arange(size_w), Numeric.Float64)
  iflag = Numeric.array([0], Numeric.Int32)
  icall = 0
  while 1:
    f = 0.
    for j in xrange(0, n, 2):
      t1 = 1.e0 - x[j]
      t2 = 1.e1 * (x[j+1] - x[j] * x[j])
      g[j+1] = 2.e1 * t2
      g[j] = -2.e0 * (x[j] * g[j+1] + t1)
      f = f + t1 * t1 + t2 * t2
    lbfgs(n, m, x, f, g, diagco, diag, iprint, eps, xtol, w, iflag)
    if (iflag[0] <= 0): break
    icall += 1
    # We allow at most 2000 evaluations of f and g
    if (icall > 2000): break

if (__name__ == "__main__"):
  import os, sys
  Endless = "--Endless" in sys.argv
  while 1:
    run()
    if (not Endless): break
  t = os.times()
  print "u+s,u,s: %.2f %.2f %.2f" % (t[0] + t[1], t[0], t[1])
