from scitbx import lbfgsb
from scitbx.array_family import flex

def exercise_minimizer_interface():
  n = 25
  m = 5
  l = flex.double(n, -1)
  u = flex.double(n, 1)
  nbd = flex.int(n)
  factr=1.0e+7
  pgtol=1.0e-5
  iprint = -1
  minimizer = lbfgsb.ext.minimizer(n, m, l, u, nbd, factr, pgtol, iprint)
  assert minimizer.task() == "START"
  x = flex.double(n, 0)
  f = 0
  g = flex.double(n, 0)
  task = minimizer.process(x, f, g)
  assert task == "FG_START"
  minimizer.request_stop()
  task = minimizer.process(x, f, g)
  assert task == "STOP: QUICK"
  minimizer.request_stop_with_restore()
  assert minimizer.task() == "STOP: CPU"

def driver1():
  n = 25
  m = 5
  nbd = flex.int(n)
  x = flex.double(n)
  l = flex.double(n)
  u = flex.double(n)
  g = flex.double(n)
  iprint = -1
  factr=1.0e+7
  pgtol=1.0e-5
  for i in xrange(0,n,2):
    nbd[i] = 2
    l[i] = 1.0e0
    u[i] = 1.0e2
  for i in xrange(1,n,2):
    nbd[i] = 2
    l[i] = -1.0e2
    u[i] = 1.0e2
  for i in xrange(n):
    x[i] = 3.0e0
  minimizer = lbfgsb.minimizer(n, m, l, u, nbd, factr, pgtol, iprint)
  f = 0
  while 0001:
    task = minimizer.process(x, f, g)
    if (task[:2] == "FG"):
      f=.25e0*(x[0]-1.e0)**2
      for i in xrange(1,n):
        f=f+(x[i]-x[i-1]**2)**2
      f=4.e0*f
      t1=x[1]-x[0]**2
      g[0]=2.e0*(x[0]-1.e0)-1.6e1*x[0]*t1
      for i in xrange(1,n-1):
        t2=t1
        t1=x[i+1]-x[i]**2
        g[i]=8.e0*t2-1.6e1*x[i]*t1
      g[n-1]=8.e0*t1
    elif (task[:5] != "NEW_X"):
      break
  assert minimizer.task() == "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

def run():
  exercise_minimizer_interface()
  driver1()
  print "OK"

if (__name__ == "__main__"):
  run()
