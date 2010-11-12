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
    exit()
  run(args=sys.argv[1:])
