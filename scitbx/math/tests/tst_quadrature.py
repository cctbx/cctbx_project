from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
import scitbx.math as sm
import math
from six.moves import range
from six.moves import zip


def twod_integrator( cub, n_points ):
  result = 0
  x = flex.double()
  y = flex.double()
  w = flex.double()
  for ii in range(n_points):
    x.append( cub.coord(ii)[0] )
    y.append( cub.coord(ii)[1] )
    w.append( cub.weight(ii) )
  tot = flex.exp( -x -y)*w
  tot = flex.sum( tot )
  return tot


def tst_cubature():
  # cubature integration of exp(-x-y)exp(-x^2-y^2)
  # 4 different cubatures are used.
  fna=sm.five_nine_1001()
  fnb=sm.five_nine_1110()
  ft=sm.seven_twelve_0120()
  nto=sm.nine_twentyone_1012()

  theory=(math.pi)*math.exp(0.5)
  r_fna=twod_integrator(fna,9)
  r_fnb=twod_integrator(fnb,9)
  r_ft=twod_integrator(ft,12)
  r_nto=twod_integrator(nto,21)
  assert approx_equal(r_fna, theory, eps=0.05)
  assert approx_equal(r_fnb, theory, eps=0.05)
  assert approx_equal(r_ft, theory, eps=0.01)
  assert approx_equal(r_nto,theory, eps=0.01)


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


def tst_gauss_legendre_engine():

  for n in range(2,95):
    gle = sm.gauss_legendre_engine(n)
    x_this = gle.x()
    step=1e-5 # check AMS pg 919. nodes for n=96 are allways more then stpe away from each other
    for ix in range( x_this.size() ):
      f = gle.f( x_this[ix] )[0]
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
  for ii in range(6,10):
    ghq = sm.gauss_hermite_engine(ii)
    w = ghq.w()
    x = ghq.x()
    f = flex.sum( (flex.exp( -x ))*w )
    assert approx_equal(f,f_theory, eps=1e-5)

  # an example of Gauss-Legendre integration
  # we integrate exp(-x) between -1 and 1
  f_theory = math.exp(1) - math.exp(-1)
  for ii in range(5,90):
    glq = sm.gauss_legendre_engine(ii)
    w = glq.w()
    x = glq.x()
    f = flex.sum( (flex.exp( -x ))*w )
    #assert approx_equal(f,f_theory, eps=1e-5)






def run():
  tst_gauss_legendre_engine()
  tst_gauss_hermite_engine()
  examples()
  tst_cubature()
  print("OK")

if (__name__ == "__main__"):
  run()
