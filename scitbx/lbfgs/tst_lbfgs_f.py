from scitbx.array_family import flex
from scitbx.lbfgs import have_lbfgs_f, lbfgs_f
from libtbx.utils import format_cpu_times
import sys

def exercise():
  n = 100
  m = 5
  x = flex.double(n)
  g = flex.double(n)
  diag = flex.double(n)
  iprint = [1, 0]
  diagco = False
  eps = 1.0e-5
  xtol = 1.0e-16
  for j in xrange(0, n, 2):
    x[j] = -1.2
    x[j+1] = 1.
  size_w = n*(2*m+1)+2*m
  w = flex.double(size_w)
  iflag = 0
  icall = 0
  while 1:
    f = 0.
    for j in xrange(0, n, 2):
      t1 = 1.e0 - x[j]
      t2 = 1.e1 * (x[j+1] - x[j] * x[j])
      g[j+1] = 2.e1 * t2
      g[j] = -2.e0 * (x[j] * g[j+1] + t1)
      f = f + t1 * t1 + t2 * t2
    iflag = lbfgs_f(
      n=n, m=m, x=x, f=f, g=g, diagco=diagco, diag=diag,
      iprint=iprint, eps=eps, xtol=xtol, w=w, iflag=iflag)
    if (iflag <= 0): break
    icall += 1
    # We allow at most 2000 evaluations of f and g
    if (icall > 2000): break

def run(args):
  if (not have_lbfgs_f):
    print "Skipping test: lbfgs_f not available."
  else:
    endless = "--endless" in args
    while 1:
      exercise()
      if (not endless): break
  print format_cpu_times()

if (__name__ == "__main__"):
  run(sys.argv[1:])
