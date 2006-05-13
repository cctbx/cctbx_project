from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
import scitbx.math as sm
import math

def tst_gauss_hermite_engine():
  # test with known values
  ghe = sm.gauss_hermite_engine(4)
  x_ams_55    = [0.524647623275290, 1.650680123885785]
  w_ams_55    = [0.8049140900055  , 0.08131283544725]
  wexs_ams_55 = [1.0599644828950  , 1.2402258176958]
  x_this = ghe.x()[0:2]
  w_this = ghe.w()[0:2]
  for x,xx in zip( x_this, x_ams_55):
    assert approx_equal( x, xx, eps=1e-8 )
  for w,ww in zip( w_this, w_ams_55):
    assert approx_equal( w, ww, eps=1e-8 )

  # test a large order set of number
  for n in range(2,29):
    ghe = sm.gauss_hermite_engine(n)
    x_this = ghe.x()
    step = 0.5/math.sqrt(n*1.0)
    for ix in range( x_this.size() ):
      f = ghe.f( x_this[ix] )[0]
      assert approx_equal(f,0,eps=1e-5)
      # check the uniqueness of each point
      for jj in  range( x_this.size() ):
        if jj != ix:
          assert ( math.fabs(x_this[ix]-x_this[jj]) >= step )

def examples():
  # an illustration of Hermite Gauss quadrature
  # we will try to integrate
  # Exp[-x^2-x] {x,-inf, inf}
  # The true answer is Exp[0.25] Sqrt[Pi]
  #
  f_theory = math.exp(0.25)*math.sqrt( math.pi )
  for ii in xrange(6,10):
    ghq = sm.gauss_hermite_engine(ii)
    w = ghq.w()
    x = ghq.x()
    f = flex.sum( (flex.exp( -x ))*w )
    assert approx_equal(f,f_theory, eps=1e-5)

def run():
  tst_gauss_hermite_engine()
  examples()
  print "OK"

if (__name__ == "__main__"):
  run()
