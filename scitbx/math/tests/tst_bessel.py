from scitbx import math
from stdlib import math as smath
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal

def exercise_interfaces():
  assert approx_equal(math.spherical_bessel(1, 2), 0.43539777498)
  assert approx_equal(math.spherical_bessel_array(1, flex.double([2])),
    [0.43539777498])

def j0(x):
  result = smath.sin(x)/x
  return result

def j1(x):
  result =  smath.sin(x)/(x*x) - smath.cos(x)/x
  return result

def j2(x):
  result = (3.0/(x*x) -1.0)*(smath.sin(x)/x) - 3.0*smath.cos(x)/(x*x)
  return result

def exercise_results():
  x = flex.double( range(1,100) )/99.0
  f0 = math.spherical_bessel_array(0,x)
  f1 = math.spherical_bessel_array(1,x)
  f2 = math.spherical_bessel_array(2,x)
  for xx, ff,fff,ffff in zip(x,f0,f1,f2):
    assert abs(ff-j0(xx))/ff < 1e-5
    assert abs(fff-j1(xx))/fff < 1e-5
    assert abs(ffff-j2(xx))/ffff < 1e-5

def tst_sph_bessel_j1():
  x = flex.double( range(1,200) )/199.0+7.5
  f1 = math.spherical_bessel_array(1,x)
  for xx,ff in zip(x,f1):
    assert abs( ff-j1(xx) )/abs(ff) < 1e-5
    #print xx, ff, j1(xx), abs( ff-j1(xx) )/abs(ff)
  print "OK"

def tst_bessel_J():
  x1 = 5.00000000e-02
  f1 = [9.99375098e-01, 2.49921883e-02, 3.12434901e-04]
  g1 = [ math.bessel_J(0,x1), math.bessel_J(1,x1), math.bessel_J(2,x1)]

  x2 = 5.00000000e-01
  f2 = [9.38469807e-01, 2.42268458e-01, 3.06040235e-02]
  g2 = [ math.bessel_J(0,x2), math.bessel_J(1,x2), math.bessel_J(2,x2) ]

  x3 = 5.00000000e+00
  f3 = [-1.77596771e-01, -3.27579138e-01, 4.65651163e-02]
  g3 = [ math.bessel_J(0,x3), math.bessel_J(1,x3), math.bessel_J(2,x3) ]

  for a,b in zip(f1,g1):
    assert (a-b)/abs(a)< 1e-5
  for a,b in zip(f2,g2):
    assert (a-b)/abs(a) < 1e-5
  for a,b in zip(f3,g3):
    assert (a-b)/abs(a) < 1e-5



def exercise():
  if (not hasattr(math, "spherical_bessel")):
    print "Skipping tst_bessel.py: functions not available."
    return
  exercise_interfaces()
  exercise_results()

def run(args):
  assert len(args) == 0
  exercise()
  print "OK"

if (__name__ == "__main__"):
  import sys
  if 'test_j1' in sys.argv[1:]:
    tst_sph_bessel_j1()
    tst_bessel_J()
    exit()
  run(args=sys.argv[1:])
