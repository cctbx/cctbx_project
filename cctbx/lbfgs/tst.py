from cctbx_boost.arraytbx import shared
from cctbx_boost import lbfgs

def run(verbose = 1):
  n = 100
  x = shared.double(n)
  g = shared.double(n)
  for j in xrange(0, n, 2):
    x[j] = -1.2
    x[j+1] = 1.
  minimizer = lbfgs.minimizer(n)
  is_converged = lbfgs.traditional_convergence_test(n)
  if (verbose):
    print "n: ", minimizer.n()
    print "m: ", minimizer.m()
    print "xtol: ", minimizer.xtol()
  while 1:
    f = 0.
    for j in xrange(0, n, 2):
      t1 = 1.e0 - x[j]
      t2 = 1.e1 * (x[j+1] - x[j] * x[j])
      g[j+1] = 2.e1 * t2
      g[j] = -2.e0 * (x[j] * g[j+1] + t1)
      f = f + t1 * t1 + t2 * t2
    if (minimizer.run(x, f, g)): continue
    if (verbose):
      print "f:", f, "gnorm:", minimizer.euclidean_norm(g)
      print minimizer.iter(), minimizer.nfun(), minimizer.stp()
    if (is_converged(x, g)): break
    if (minimizer.nfun() > 2000): break
    assert minimizer.run(x, f, g)

if (__name__ == "__main__"):
  import os, sys
  Endless = "--Endless" in sys.argv
  while 1:
    run(not Endless)
    if (not Endless): break
  t = os.times()
  print "u+s,u,s:", t[0] + t[1], t[0], t[1]
