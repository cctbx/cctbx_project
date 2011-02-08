import scitbx.math
import boost.rational
from scitbx.math import line_given_points
from scitbx.math import dihedral_angle
from scitbx.math import euler_angles_as_matrix
from scitbx.math import erf_verification, erf, erfc, erfcx
from scitbx.math import bessel_i1_over_i0, bessel_i0, bessel_i1,\
     bessel_ln_of_i0, ei1, ei0
from scitbx.math import bessel_inverse_i1_over_i0
from scitbx.math import gamma_incomplete, gamma_incomplete_complement
from scitbx.math import gamma_complete, exponential_integral_e1z
from scitbx.math import lambertw
from scitbx.math import golay_24_12_generator
from scitbx.math import inertia_tensor
from scitbx.math import principal_axes_of_inertia, \
     principal_axes_of_inertia_2d
from scitbx.math import sphere_3d, minimum_covering_sphere
from scitbx.math import signed_phase_error, phase_error, nearest_phase
from scitbx.math import icosahedron
from scitbx.math import chebyshev_base
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_fitter
from scitbx.math import slatec_dgamma, slatec_dlngam
from scitbx.math import distributions
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.utils import user_plus_sys_time
from libtbx.test_utils import Exception_expected, approx_equal, eps_eq
import libtbx.load_env
from cStringIO import StringIO
from itertools import count
import random
import math
import time
import sys

if (libtbx.env.has_module("tntbx")):
  import tntbx
else:
  tntbx = None

def exercise_div_mod():
  from scitbx.math import divmod
  q,r = divmod(4.34, 2)
  assert q == 2
  assert approx_equal(r, 0.34)
  q,r = divmod(-3.23, 2)
  assert q == -2
  assert approx_equal(r, 0.77)

  for y,q,r in [(5.23, 2, 1.3), (-5.23, 2, 1.3),
                (5.23, -2, 1.3), (-5.23, -2, 1.3)]:
    x = q*y + r
    q1,r1 = divmod(x,y)
    assert q1 == q
    assert approx_equal(r1,r)

  assert divmod(0, 1) == (0,0)
  assert divmod(1, 1) == (1,0)
  assert divmod(2, 1) == (2,0)
  assert divmod(-0, 1) == (0,-0)
  assert divmod(-1, 1) == (-1,-0)
  assert divmod(-2, 1) == (-2,-0)

  assert divmod(0, -1) == (0,0)
  assert divmod(1, -1) == (-1,0)
  assert divmod(2, -1) == (-2,0)
  assert divmod(-0, -1) == (0,-0)
  assert divmod(-1, -1) == (1,-0)
  assert divmod(-2, -1) == (2,-0)

  assert divmod(1, 2) == (0,1)
  assert divmod(3, 2) == (2,-1)
  assert divmod(5, 2) == (2,1)
  assert divmod(-1, 2) == (0,-1)
  assert divmod(-3, 2) == (-2,1)
  assert divmod(-5, 2) == (-2,-1)

  assert divmod(1, -2) == (0,1)
  assert divmod(3, -2) == (-2,-1)
  assert divmod(5, -2) == (-2,1)
  assert divmod(-1, -2) == (0,-1)
  assert divmod(-3, -2) == (2,1)
  assert divmod(-5, -2) == (2,-1)

def exercise_floating_point_epsilon():
  float_eps = scitbx.math.floating_point_epsilon_float_get()
  double_eps = scitbx.math.floating_point_epsilon_double_get()
  assert 1.+float_eps != 1.
  assert 1.+double_eps != 1.
  assert float_eps >= double_eps
  assert 1.+double_eps/2. == 1.

def exercise_line_given_points():
  lgp = line_given_points(points=[(0,0,0), (0,0,0)])
  assert lgp.distance_sq(point=matrix.col((0,0,0))) == 0
  assert lgp.distance_sq(point=matrix.col((0,0,1))) == 1
  lgp = line_given_points(points=[(0,0,1), (0,0,0)])
  assert lgp.distance_sq(point=matrix.col((0,0,1))) == 0
  assert lgp.distance_sq(point=matrix.col((0,1,0))) == 1
  lgp = line_given_points(points=[(1,2,3), (3,1,4)])
  assert lgp.distance_sq(point=matrix.col((0,0,0))) == 12
  assert lgp.distance_sq(point=matrix.col((0,0,1))) == 8

def exercise_dihedral_angle():
  def dihe(sites):
    t = dihedral_angle(sites=sites, deg=True)
    v = dihedral_angle(sites=flex.vec3_double(sites), deg=True)
    if (t is None):
      assert v is None
    else:
      assert v is not None
      assert approx_equal(t, v)
    return t
  d = dihe(sites=[(0,0,0)]*4)
  assert d is None
  d = dihe(sites=[(1,0,0), (0,0,0), (0,1,0), (1,1,0)])
  assert approx_equal(d, 0)
  d = dihe(sites=[(1,0,0), (0,0,0), (0,1,0), (-1,1,0)])
  assert approx_equal(abs(d), 180)
  d = dihe(sites=[(1,0,0), (0,0,0), (0,1,0), (0,1,1)])
  assert approx_equal(d, -90)
  d = dihe(sites=[(1,0,0), (0,0,0), (0,1,0), (0,1,-1)])
  assert approx_equal(d, 90)
  mt = flex.mersenne_twister(seed=0)
  for ad in xrange(-179, 180):
    ar = math.radians(ad)
    c, s = math.cos(ar), math.sin(ar)
    f = 2 - mt.random_double()
    sites = flex.vec3_double([(f,0,0), (0,0,0), (0,f,0), (f*c,f,-f*s)])
      #                                                         ^
      #                                 conventions are wonderful
    d = dihedral_angle(sites=sites, deg=True)
    assert approx_equal(d, ad)
    r = mt.random_double_r3_rotation_matrix()
    t = mt.random_double_point_on_sphere()
    d = dihedral_angle(sites=r*sites+t, deg=True)
    assert approx_equal(d, ad)
    d = matrix._dihedral_angle(sites=matrix.col_list(r*sites+t), deg=True)
    assert approx_equal(d, ad)
  sites = [
    (-3.193, 1.904, 4.589),
    (-1.955, 1.332, 3.895),
    (-1.005, 2.228, 3.598),
    ( 0.384, 1.888, 3.199)]
  d = dihe(sites=sites)
  assert approx_equal(d, 166.212120415)

def exercise_euler_angles():
  assert approx_equal(euler_angles_as_matrix([0,0,0]).elems,
    [1,0,0,0,1,0,0,0,1])
  assert approx_equal(euler_angles_as_matrix([180,0,0], deg=True).elems,
    [-1,0,0,0,-1,0,0,0,1])
  assert approx_equal(euler_angles_as_matrix([0,180,0], deg=True).elems,
    [-1,0,0,0,1,0,0,0,-1])
  assert approx_equal(euler_angles_as_matrix([0,0,180], deg=True).elems,
    [-1,0,0,0,-1,0,0,0,1])
  assert approx_equal(euler_angles_as_matrix([90,0,0], deg=True).elems,
    [0,-1,0,1,0,0,0,0,1])
  assert approx_equal(euler_angles_as_matrix([0,90,0], deg=True).elems,
    [0,0,1,0,1,0,-1,0,0])
  assert approx_equal(euler_angles_as_matrix([0,0,90], deg=True).elems,
    [0,-1,0,1,0,0,0,0,1])

def exercise_erf():
  erf_verify = erf_verification()
  # Expected results obtained by running the original
  # FORTRAN code under Tru64 Unix.
  # See also: scitbx/math/dev/README
  erf_verify(erf, 1.102483591684125E-004, 1.244019511880774E-004)
  erf_verify(erf, 0.117651583521854, 0.132145601355679)
  erf_verify(erf, 0.235069688589948, 0.260442011292979)
  erf_verify(erf, 0.352309859191852, 0.381686081076891)
  erf_verify(erfc, 0.970883582602861, 0.169740928450908)
  erf_verify(erfc, 0.662251325736969, 0.348982462199060)
  erf_verify(erfc, 1.35402517159916, 5.550771183345651E-002)
  erf_verify(erfc, 1.04576304738078, 0.139158413714973)
  erf_verify(erfc, 1.73790185144654, 1.398048693929801E-002)
  erf_verify(erfc, 1.42912826254883, 4.327018321700024E-002)
  erf_verify(erfc, 2.12125887397317, 2.700566694860312E-003)
  erf_verify(erfc, 1.81249050086459, 1.036977546108103E-002)
  erf_verify(erfcx, 0.973392426408860, 0.434963885763398)
  erf_verify(erfcx, 0.665081793626574, 0.539933818979981)
  erf_verify(erfcx, 1.35724600005572, 0.346596780023340)
  erf_verify(erfcx, 1.04836345595439, 0.414722416503421)
  erf_verify(erfcx, 1.74107610470511, 0.286145510591003)
  erf_verify(erfcx, 1.43212197979175, 0.333053266914237)
  erf_verify(erfcx, 2.12428846184817, 0.242739423717836)
  erf_verify(erfcx, 1.81606172854839, 0.276556810628268)
  erf_verify(erfc, 2.62003792003157, 2.111463761469150E-004)
  erf_verify(erfc, 5.12059159283289, 4.433890735237740E-013)
  erf_verify(erfc, 8.63668354446932, 2.613296860860847E-034)
  erf_verify(erfc, 11.1437916191963, 5.891081160010891E-056)
  erf_verify(erfc, 14.6497556694182, 2.389628051234427E-095)
  erf_verify(erfc, 17.1517732073587, 5.677885934529786E-130)
  erf_verify(erfc, 20.6668792553407, 8.706194441633715E-188)
  erf_verify(erfc, 23.1638554546886, 2.287385773272363E-235)
  erf_verify(erfcx, 2.62907657906194, 0.201608800911889)
  erf_verify(erfcx, 4.37646062474290, 0.125783546292562)
  erf_verify(erfcx, 7.13820945795942, 7.828417109269534E-002)
  erf_verify(erfcx, 8.88772516645629, 6.308522532999319E-002)
  erf_verify(erfcx, 11.6501072021873, 4.825137678138556E-002)
  erf_verify(erfcx, 13.4008106775670, 4.198489891639072E-002)
  erf_verify(erfcx, 16.1586604627164, 3.484913452473175E-002)
  erf_verify(erfcx, 17.9110272449630, 3.145069912820161E-002)
  erf_verify(erfcx, 20.4992639479713, 2.748979983772347E-002)
  erf_verify(erfcx, 19.9992639479713, 2.817538309213394E-002)
  erf_verify(erf, 0.000000000000000E+000, 0.000000000000000E+000)
  erf_verify(erf, 0.000000000000000E+000, 0.000000000000000E+000)
  erf_verify(erfc, 0.000000000000000E+000, 1.00000000000000)
  erf_verify(erfcx, 0.000000000000000E+000, 1.00000000000000)
  erf_verify(erf, -0.500000000000000, -0.520499877813047)
  erf_verify(erf, 0.500000000000000, 0.520499877813047)
  erf_verify(erfc, -0.500000000000000, 1.52049987781305)
  erf_verify(erfcx, -0.500000000000000, 1.95236048918256)
  erf_verify(erf, -1.00000000000000, -0.842700792949715)
  erf_verify(erf, 1.00000000000000, 0.842700792949715)
  erf_verify(erfc, -1.00000000000000, 1.84270079294971)
  erf_verify(erfcx, -1.00000000000000, 5.00898008076228)
  erf_verify(erf, -1.50000000000000, -0.966105146475311)
  erf_verify(erf, 1.50000000000000, 0.966105146475311)
  erf_verify(erfc, -1.50000000000000, 1.96610514647531)
  erf_verify(erfcx, -1.50000000000000, 18.6538862562627)
  erf_verify(erf, -2.00000000000000, -0.995322265018953)
  erf_verify(erf, 2.00000000000000, 0.995322265018953)
  erf_verify(erfc, -2.00000000000000, 1.99532226501895)
  erf_verify(erfcx, -2.00000000000000, 108.940904389978)
  erf_verify(erf, -2.50000000000000, -0.999593047982555)
  erf_verify(erf, 2.50000000000000, 0.999593047982555)
  erf_verify(erfc, -2.50000000000000, 1.99959304798255)
  erf_verify(erfcx, -2.50000000000000, 1035.81484297262)
  erf_verify(erf, -3.00000000000000, -0.999977909503001)
  erf_verify(erf, 3.00000000000000, 0.999977909503001)
  erf_verify(erfc, -3.00000000000000, 1.99997790950300)
  erf_verify(erfcx, -3.00000000000000, 16205.9888539996)
  erf_verify(erf, -3.50000000000000, -0.999999256901628)
  erf_verify(erf, 3.50000000000000, 0.999999256901628)
  erf_verify(erfc, -3.50000000000000, 1.99999925690163)
  erf_verify(erfcx, -3.50000000000000, 417962.422445770)
  erf_verify(erf, -4.00000000000000, -0.999999984582742)
  erf_verify(erf, 4.00000000000000, 0.999999984582742)
  erf_verify(erfc, -4.00000000000000, 1.99999998458274)
  erf_verify(erfcx, -4.00000000000000, 17772220.9040163)
  erf_verify(erf, -4.50000000000000, -0.999999999803384)
  erf_verify(erf, 4.50000000000000, 0.999999999803384)
  erf_verify(erfc, -4.50000000000000, 1.99999999980338)
  erf_verify(erfcx, -4.50000000000000, 1245928884.27441)
  erf_verify(erf, 1.797693134862316E+308, 1.00000000000000)
  erf_verify(erf, 0.000000000000000E+000, 0.000000000000000E+000)
  erf_verify(erfc, 0.000000000000000E+000, 1.00000000000000)
  erf_verify(erfc, -1.797693134862316E+308, 2.00000000000000)
  erf_verify(erfc, 19.9074438229742, 2.178791258493533E-174)
  erf_verify(erfc, 26.5432584306324, 0.000000000000000E+000)
  erf_verify(erfcx, 2.377124393213978E+307, 2.373412115741015E-308)
  erf_verify(erfcx, -23.9658621423763, 5.540070644707187E+249)
  erf_verify(erfcx, -26.6287357137515, 1.790000000000000E+308)
  assert erf_verify.max_delta < erf_verify.tolerance

def exercise_exponential_integral_e1z():
  assert approx_equal(exponential_integral_e1z(0.5), 0.559773595)
  assert approx_equal(exponential_integral_e1z(1.0), 0.219383934)
  assert approx_equal(exponential_integral_e1z(1.5), 0.100019582)
  assert approx_equal(exponential_integral_e1z(2.0), 0.048900511)


def exercise_gamma_incomplete():
  assert approx_equal(gamma_incomplete(2.0, 0.1),0.004678840160445)
  assert approx_equal(gamma_incomplete(2.0, 0.5),0.090204010431050)
  assert approx_equal(gamma_incomplete(2.0, 2.5),0.712702504816354)
  assert approx_equal(gamma_incomplete(2.0, 5.0),0.959572318005487)
  assert approx_equal(gamma_incomplete(2.0,15.5),0.999996938604252)
  assert approx_equal(gamma_incomplete(2.0,21.0),0.999999983318367)
  #
  assert approx_equal(gamma_incomplete(20.0, 0.1),0)
  assert approx_equal(gamma_incomplete(20.0, 0.5),0)
  assert eps_eq(gamma_incomplete(20.0, 2.5),
    3.480438159897403e-12,eps=1.e-5)
  assert eps_eq(gamma_incomplete(20.0, 5.0),3.452135821646607e-7)
  assert approx_equal(gamma_incomplete(20.0,15.5),0.154492096867129)
  assert approx_equal(gamma_incomplete(20.0,21.0),0.615737227735658)
  try: gamma_incomplete(a=20.0, x=15.5, max_iterations=5)
  except RuntimeError, e:
    assert str(e) == \
      "scitbx Error: gamma::incomplete_series(" \
      "a=20, x=15.5, max_iterations=5) failed to converge"
  else: raise Exception_expected
  try: gamma_incomplete(a=20.0, x=25.5, max_iterations=5)
  except RuntimeError, e:
    assert str(e) == \
      "scitbx Error: gamma::incomplete_continued_fraction(" \
      "a=20, x=25.5, max_iterations=5) failed to converge"
  else: raise Exception_expected
  #
  assert approx_equal(1-gamma_incomplete_complement(2.0, 2.5),
                                   gamma_incomplete(2.0, 2.5))

def exercise_gamma_complete():
  ## complete gamma with lanczos approx for x<12 and minimax otherwise
  assert approx_equal(gamma_complete(0.1),9.5135076986687)
  assert approx_equal(gamma_complete(0.5),1.7724538509055)
  assert approx_equal(gamma_complete(2.5),1.3293403881791)
  assert approx_equal(gamma_complete(5.0),24.0)
  assert approx_equal(gamma_complete(15.5),3.3483860987356E11)
  assert approx_equal(gamma_complete(21.0),2432902008176640000)
  ## complete gamma with lanczos approx for all values
  assert approx_equal(gamma_complete(0.1,minimax=False),9.5135076986687)
  assert approx_equal(gamma_complete(0.5,minimax=False),1.7724538509055)
  assert approx_equal(gamma_complete(2.5,minimax=False),1.3293403881791)
  assert approx_equal(gamma_complete(5.0,minimax=False),24.)
  assert approx_equal(gamma_complete(15.5,minimax=False),3.3483860987356E11)
  assert approx_equal(gamma_complete(21.0,minimax=False),2432902008176640000)
  assert "%.8g" % gamma_complete(171.624-1.e-6) == "1.7942025e+308"
  #
  try: gamma_complete(171.624)
  except RuntimeError, e:
    assert str(e) \
        == "scitbx Error: gamma::complete_minimax(171.624): domain error"
  else: raise Exception_expected
  assert "%.8g" % gamma_complete(141.691-1.e-6) == "4.1104518e+242"
  try: gamma_complete(141.691, minimax=False)
  except RuntimeError, e:
    assert str(e) \
        == "scitbx Error: gamma::complete_lanczos(141.691): domain error"
  else: raise Exception_expected

def exercise_bessel():
  assert approx_equal(bessel_i1_over_i0(-1e+9), -1.0)
  assert approx_equal(bessel_i1_over_i0(-99.99),-0.994988)
  assert approx_equal(bessel_i1_over_i0(-50.00),-0.98995)
  assert approx_equal(bessel_i1_over_i0( -1.00),-0.44639)
  assert approx_equal(bessel_i1_over_i0(  0.0),  0.0)
  assert approx_equal(bessel_i1_over_i0(  1.00), 0.44639)
  assert approx_equal(bessel_i1_over_i0( 50.0),  0.98995)
  assert approx_equal(bessel_i1_over_i0( 99.99), 0.994988)
  assert approx_equal(bessel_i1_over_i0(1e+9),   1.0)
  x=0.0
  while x <= 100.0:
    assert approx_equal(-bessel_i1_over_i0(-x),bessel_i1_over_i0(x))
    x+=0.01
  x=0.0
  while x <= 10.0:
    assert eps_eq(x,bessel_inverse_i1_over_i0(bessel_i1_over_i0(x)),eps=5.e-2)
    x+=0.01
  assert approx_equal(bessel_i0(0.0), 1.0)
  assert approx_equal(bessel_i1(0.0), 0.0)
  x=-500.0
  while x <= 500.0:
    a = bessel_i1_over_i0(x)
    b = bessel_i1(x) / bessel_i0(x)
    assert approx_equal(a,b,1.e-5)
    x+=0.01
  x=-100.0
  while x <= 100.0:
    assert approx_equal(bessel_ln_of_i0(x),math.log(bessel_i0(x)))
    x+=0.01

def exercise_eix():
  x = flex.double( range(4000) )/20.0
  expx = flex.exp( -x )
  for xx, ex in zip(x,expx):
    tmp_i0  = bessel_i0(xx)
    tmp_i1  = bessel_i1(xx)
    tmp_ei0 = ei0(xx)
    tmp_ei1 = ei1(xx)
    assert approx_equal(  tmp_i0*ex, tmp_ei0, eps=1e-3 )
    assert approx_equal(  tmp_i1*ex, tmp_ei1, eps=1e-3 )

def exercise_random_cheb_polynome(n_terms,
                             low_limit,
                             high_limit,
                             h=0.00001):
  x = flex.double(range(100))/101.0
  x = x*(high_limit-low_limit) + low_limit
  coefs = flex.random_double(n_terms)
  cheb = chebyshev_polynome(n_terms,
                            low_limit,
                            high_limit,
                            coefs)
  y = cheb.f(x)
  y_tmp = cheb.f(x+h)
  dydx = cheb.dfdx(x)
  for ii in range(100):
    assert approx_equal(dydx[ii], (y_tmp[ii]-y[ii])/h,eps=1e-4)


def exercise_cheb_fitter(n_terms,
                    low_limit,
                    high_limit,
                    h=0.0000001):
  x = flex.double(range(100))/101.0
  x = x*(high_limit-low_limit) + low_limit
  coefs = flex.random_double(n_terms)
  cheb_fitter = chebyshev_fitter(n_terms,
                                 low_limit,
                                 high_limit,
                                 coefs)

  finite_diffs = flex.double(n_terms,0)

  for ii in range(100):
    exact = cheb_fitter.dfdcoefs(x[ii])
    f = cheb_fitter.f(x[ii])
    for jj in range(n_terms):
      coefs[jj]+=h
      cheb_fitter.replace(coefs)
      df = cheb_fitter.f(x[ii])
      coefs[jj]-=h
      finite_diffs[jj] = (df-f)/h
      cheb_fitter.replace(coefs)
      assert approx_equal(exact[jj], finite_diffs[jj])


def exercise_cheb_base_and_polynome():
  x = flex.double(100,0)
  y1 = flex.double(100,0)
  y2 = flex.double(100,0)
  y3 = flex.double(100,0)
  dy1 = flex.double(100,0)
  dy2 = flex.double(100,0)
  dy3 = flex.double(100,0)


  cheb_1_coefs = flex.double([0,1])
  cheb_1 = chebyshev_base(2,-1.,1.,cheb_1_coefs )
  cheb_1_d = chebyshev_polynome(2,-1.,1.,cheb_1_coefs )

  cheb_2_coefs = flex.double([0,0,1])
  cheb_2 = chebyshev_base(3,-1.,1.,cheb_2_coefs )
  cheb_2_d = chebyshev_polynome(3,-1.,1.,cheb_2_coefs )

  cheb_3_coefs = flex.double([0,0,0,1])
  cheb_3 = chebyshev_base(4,-1.,1.,cheb_3_coefs )
  cheb_3_d = chebyshev_polynome(4,-1.,1.,cheb_3_coefs )


  for ii in range(100):
    x[ii] = (ii-50)/51.0

    y1[ii] = x[ii]
    dy1[ii]= 1.0

    y2[ii] = 2.0* x[ii]*x[ii] - 1.0
    dy2[ii]=4.0*x[ii]

    y3[ii] = 4.0*x[ii]*x[ii]*x[ii]-3.0*x[ii]
    dy3[ii] = 12.0*x[ii]*x[ii] - 3.0


  cheb_1_y = cheb_1.f(x)
  cheb_2_y = cheb_2.f(x)
  cheb_3_y = cheb_3.f(x)

  cheb_1_d_y = cheb_1_d.f(x)
  cheb_2_d_y = cheb_2_d.f(x)
  cheb_3_d_y = cheb_3_d.f(x)

  cheb_1_d_ydx = cheb_1_d.dfdx(x)
  cheb_2_d_ydx = cheb_2_d.dfdx(x)
  cheb_3_d_ydx = cheb_3_d.dfdx(x)

  for ii in range(100):
    assert approx_equal(y1[ii],cheb_1_y[ii],eps=1e-4)
    assert approx_equal(y1[ii],cheb_1_d_y[ii],eps=1e-4)
    assert approx_equal(dy1[ii],cheb_1_d_ydx[ii],eps=1e-4)

    assert approx_equal(y2[ii],cheb_2_y[ii],eps=1e-4)
    assert approx_equal(y2[ii],cheb_2_d_y[ii],eps=1e-4)
    assert approx_equal(dy2[ii],cheb_2_d_ydx[ii],eps=1e-4)

    assert approx_equal(y2[ii],cheb_2_y[ii],eps=1e-4)
    assert approx_equal(y2[ii],cheb_2_d_y[ii],eps=1e-4)
    assert approx_equal(dy2[ii],cheb_2_d_ydx[ii],eps=1e-4)

def exercise_cheb_family():
  exercise_cheb_base_and_polynome()
  for ii in range(10):
    exercise_cheb_fitter(5,-1,1,0.00000001)
    exercise_cheb_fitter(5,-10,1,0.00000001)
    exercise_cheb_fitter(5,-10,10,0.00000001)
    exercise_cheb_fitter(8,-1,1,0.00000001)
    exercise_cheb_fitter(8,-10,1,0.00000001)
    exercise_cheb_fitter(8,-10,10,0.00000001)
    exercise_random_cheb_polynome(8, -1.0, 1.0,0.000000001)
    exercise_random_cheb_polynome(5, -1.0, 1.0,0.000000001)
    exercise_random_cheb_polynome(8, 1.0, 10.0,0.000000001)
    exercise_random_cheb_polynome(5, -10.0, 1.0,0.000000001)

def check_lambertw(x):
  w = lambertw(x=x)
  assert eps_eq(w*math.exp(w), x)

def exercise_lambertw():
  check_lambertw(-math.exp(-1)+1.e-4)
  check_lambertw(-1.e-5)
  check_lambertw(0)
  check_lambertw(1.e-5)
  check_lambertw(1-1.e-5)
  check_lambertw(1+1.e-5)
  check_lambertw(3-1.e-5)
  check_lambertw(3+1.e-5)
  for i in xrange(100):
    check_lambertw(x=i/10.-0.35)
  for i in xrange(20):
    check_lambertw(x=2.**i)
    check_lambertw(x=5.**i)
    check_lambertw(x=10.**i)
  try: lambertw(x=-math.exp(-1)-1.e-4)
  except RuntimeError, e:
    assert str(e) == "lambertw(x) domain error: x < -exp(-1)"
  else: raise Exception_expected
  try: lambertw(x=1, max_iterations=1)
  except RuntimeError, e:
    assert str(e) == "lambertw error: iteration did not converge"
  else: raise Exception_expected

def matrix_mul(a, ar, ac, b, br, bc):
  assert br == ac
  result = []
  for i in xrange(ar):
    for k in xrange(bc):
      s = 0
      for j in xrange(ac):
        s += a[i * ac + j] * b[j * bc + k]
      result.append(s)
  return result

def exercise_golay():
  weights = [0]*25
  gg = golay_24_12_generator()
  while not gg.at_end():
    weights[list(gg.next()).count(1)] += 1
  assert weights == [1,0,0,0,0,0,0,0,759,0,0,0,2576,0,0,0,759,0,0,0,0,0,0,0,1]
  try:
    gg.next()
  except StopIteration, e:
    assert str(e) == "golay_24_12_generator is exhausted."
  else:
    raise Exception_expected
  weights = [0]*25
  for code in golay_24_12_generator():
    weights[list(code).count(1)] += 1
  assert weights == [1,0,0,0,0,0,0,0,759,0,0,0,2576,0,0,0,759,0,0,0,0,0,0,0,1]

def exercise_inertia_tensor():
  t = inertia_tensor(
    points=flex.vec3_double(), pivot=(0,0,0))
  assert t == (0,0,0,0,0,0)
  t = inertia_tensor(
    points=flex.vec3_double(), weights=flex.double(), pivot=(0,0,0))
  assert t == (0,0,0,0,0,0)
  t = inertia_tensor(
    points=flex.vec3_double([(4,6,2),(2,6,4)]), pivot=(-1,8,4))
  assert approx_equal(t, [12, 38, 42, 16, 10, -4])
  t = inertia_tensor(
    points=flex.vec3_double([(4,6,2),(2,6,4)]),
    weights=flex.double([2,3]),
    pivot=(-1,8,4))
  assert approx_equal(t, [28, 85, 97, 38, 20, -8])

def exercise_principal_axes_of_inertia():
  rnd = random.random
  for i_trial in xrange(10):
    if (i_trial == 0):
      points = flex.vec3_double()
    elif (i_trial == 1):
      points = flex.vec3_double([[0,0,0]])
    else:
      points = flex.vec3_double([[rnd(),rnd(),rnd()]])
    pai = principal_axes_of_inertia(points=points)
    if (i_trial == 0):
      assert approx_equal(pai.center_of_mass(), [0,0,0])
    else:
      assert approx_equal(pai.center_of_mass(), points[0])
    assert approx_equal(pai.inertia_tensor(), [0,0,0,0,0,0])
    es = pai.eigensystem()
    assert approx_equal(es.values(), [0,0,0])
    assert approx_equal(es.vectors(), [1,0,0,0,1,0,0,0,1])
    assert approx_equal(
      pai.change_of_basis_mx_to_principal(), [1,0,0,0,1,0,0,0,1])
    assert pai.distance_to_inertia_ellipsoid_surface(
      unit_direction=(1,0,0)) == 0
  for i_trial in xrange(10):
    if (i_trial == 0):
      center_of_mass = [0,0,0]
    else:
      center_of_mass = [rnd(),rnd(),rnd()]
    points = flex.vec3_double(flex.nested_loop([-1,-1,-1], [2,2,2]))
    points += center_of_mass
    pai = principal_axes_of_inertia(points=points)
    assert approx_equal(pai.center_of_mass(), center_of_mass)
    assert approx_equal(pai.inertia_tensor(), [36,36,36,0,0,0])
    es = pai.eigensystem()
    assert approx_equal(es.values(), [36,36,36])
    cp = pai.change_of_basis_mx_to_principal()
    assert approx_equal(matrix.sqr(cp).determinant(), 1)
    vectors = [matrix.col(es.vectors()[i:i+3]) for i in range(3)]
    # testing for specific eigenvectors for the case
    # of degenerate eigenvalues is a bit risky
    paip = principal_axes_of_inertia(points=cp*points)
    assert approx_equal(paip.inertia_tensor(), [36,36,36,0,0,0])
    assert approx_equal(pai.distance_to_inertia_ellipsoid_surface(
      unit_direction=(1,0,0)), 36)
    assert pai.distance_to_inertia_ellipsoid_surface(
      unit_direction=(0,0,0)) == 0
  for i_trial in xrange(10):
    # test for the case of non-degenerate eigenvalues
    # check that the inertia tensor and eigenvectors
    # transform correctly under rotation
    eps = 1e-12
    if (i_trial == 0):
      center_of_mass = [0,0,0]
      rotation = (1,0,0,0,1,0,0,0,1)
    else:
      center_of_mass = [rnd(),rnd(),rnd()]
      rotation = flex.random_double_r3_rotation_matrix()
    # a parallelepiped
    points = flex.vec3_double([
      (-4,-2,-1), (-3,-2,-1), (-2,-2,-1),
      (-3,-1,-1), (-2,-1,-1), (-1,-1,-1),
      (-2, 0,-1), (-1, 0,-1), ( 0, 0,-1),
      (-2,-1, 0), (-1,-1, 0), ( 0,-1, 0),
      (-1, 0, 0), ( 0, 0, 0), ( 1, 0, 0),
      ( 0, 1, 0), ( 1, 1, 0), ( 2, 1, 0),
      ( 0, 0, 1), ( 1, 0, 1), ( 2, 0, 1),
      ( 1, 1, 1), ( 2, 1, 1), ( 3, 1, 1),
      ( 2, 2, 1), ( 3, 2, 1), ( 4, 2, 1),
      ])
    points = rotation * points
    points += center_of_mass
    pai = principal_axes_of_inertia(points=points)
    es = pai.eigensystem()
    assert approx_equal(pai.center_of_mass(), center_of_mass)
    R = matrix.sqr(rotation)
    R_t = R.transpose()
    assert approx_equal(
      (R_t * matrix.sym(sym_mat3=pai.inertia_tensor()) * R).as_sym_mat3(),
      (54,126,144,-54,-36,-18), eps=eps)
    assert approx_equal(es.values(),
      [156.90386550855695, 154.33160314031173, 12.764531351131396], eps=eps)
    expected_vectors = [
      (-0.44909878511104717, 0.29312841385740002, 0.84402962874606857),
      (-0.29312841385727212, 0.84402962874598531, -0.44909878511128715),
      (0.84402962874611298, 0.4490987851112036, 0.29312841385703242)]
    for i in range(3):
      vec = matrix.col(es.vectors()[3*i:3*(i+1)])
      try:
        assert approx_equal(R_t * vec, expected_vectors[i], eps=eps, out=None)
      except AssertionError:
        # we don't know the direction of the eigenvector
        assert approx_equal(- R_t * vec, expected_vectors[i], eps=eps)
  for i_trial in xrange(10):
    if (i_trial == 0):
      center_of_mass = [0,0,0]
    else:
      center_of_mass = [rnd(),rnd(),rnd()]
    points = flex.vec3_double()
    for point in flex.nested_loop([-1,-1,-1], [2,2,2]):
      points.append((matrix.col([point[0],point[1]*2,point[2]*3])
                   + matrix.col(center_of_mass)).elems)
    pai = principal_axes_of_inertia(points=points)
    assert approx_equal(pai.center_of_mass(), center_of_mass)
    assert approx_equal(pai.inertia_tensor(), [234,180,90,0,0,0])
    es = pai.eigensystem()
    assert approx_equal(es.values(), [234,180,90])
    if (i_trial == 0):
      assert approx_equal(es.vectors(), [1,0,0,0,1,0,0,0,1])
    else:
      cp = pai.change_of_basis_mx_to_principal()
      assert approx_equal(matrix.sqr(cp).determinant(), 1)
      paip = principal_axes_of_inertia(points=cp*points)
      assert approx_equal(paip.inertia_tensor(), [234,180,90,0,0,0])
  for i_trial in xrange(10):
    if (i_trial == 0):
      center_of_mass = [0,0,0]
    else:
      center_of_mass = [rnd(),rnd(),rnd()]
    if (i_trial < 2):
      rot = matrix.sqr([1,0,0,0,1,0,0,0,1])
    else:
      rot = euler_angles_as_matrix(
        angles=[random.uniform(0,360) for i in xrange(3)],
        deg=True)
    points = flex.vec3_double()
    for point in [
      [-1,-1, 0],[-1, 1, 0],
      [ 1,-1, 0],[ 1, 1, 0],
      [-1, 0,-1],[-1, 0, 0],[-1, 0, 1],
      [ 1, 0,-1],[ 1, 0, 0],[ 1, 0, 1],
      [ 0, 0,-1],
      [ 0, 0, 1]]:
      points.append((rot*matrix.col(point)+matrix.col(center_of_mass)).elems)
    pai = principal_axes_of_inertia(points=points)
    assert approx_equal(pai.center_of_mass(), center_of_mass)
    if (i_trial < 2):
      assert approx_equal(pai.inertia_tensor(), [10,16,14,0,0,0])
    es = pai.eigensystem()
    assert approx_equal(es.values(), [16,14,10])
    if (i_trial < 2):
      assert approx_equal(es.vectors(), [0,1,0,0,0,1,1,0,0])
    else:
      cp = pai.change_of_basis_mx_to_principal()
      assert approx_equal(matrix.sqr(cp).determinant(), 1)
      paip = principal_axes_of_inertia(points=cp*points)
      assert approx_equal(paip.inertia_tensor(), [16,14,10,0,0,0])
    assert abs(abs(matrix.col(es.vectors()[0:3]).dot(
                   rot*matrix.col([0,1,0])))-1) < 1.e-3
    assert abs(abs(matrix.col(es.vectors()[3:6]).dot(
                   rot*matrix.col([0,0,1])))-1) < 1.e-3
    assert abs(abs(matrix.col(es.vectors()[6:9]).dot(
                   rot*matrix.col([1,0,0])))-1) < 1.e-3
    weights = flex.double(points.size(), 1)
    points_plus = points.deep_copy()
    for i_p in [0,3,5,5,9,10,10,10]:
      points_plus.append(points[i_p])
      weights[i_p] += 1
    paip = principal_axes_of_inertia(points=points_plus)
    paiw = principal_axes_of_inertia(points=points, weights=weights)
    assert approx_equal(paip.center_of_mass(), paiw.center_of_mass())
    assert approx_equal(paip.inertia_tensor(), paiw.inertia_tensor())

def exercise_principal_axes_of_inertia_2d():
  rnd = random.random
  for i_trial in xrange(10):
    if (i_trial == 0):
      points = flex.vec2_double()
    elif (i_trial == 1):
      points = flex.vec2_double([[0,0]])
    else:
      points = flex.vec2_double([[rnd(),rnd()]])
    pai = principal_axes_of_inertia_2d(points=points)
    if (i_trial == 0):
      assert approx_equal(pai.center_of_mass(), [0,0])
    else:
      assert approx_equal(pai.center_of_mass(), points[0])
    assert approx_equal(pai.inertia_tensor(), [0,0,0])
    es = pai.eigensystem()
    assert approx_equal(es.values(), [0,0])
    assert approx_equal(es.vectors(), [1,0,0,1])
    assert pai.distance_to_inertia_ellipsoid_surface(
      unit_direction=(1,0)) == 0
  for i_trial in xrange(10):
    if (i_trial == 0):
      center_of_mass = [0,0]
    else:
      center_of_mass = [rnd(),rnd()]
    points = flex.vec2_double()
    for point in flex.nested_loop([-1,-1], [2,2]):
      points.append((matrix.col(point) + matrix.col(center_of_mass)).elems)
    pai = principal_axes_of_inertia_2d(points=points)
    assert approx_equal(pai.center_of_mass(), center_of_mass)
    assert approx_equal(pai.inertia_tensor(), [6,6,0])
    es = pai.eigensystem()
    assert approx_equal(es.values(), [6,6])
    if (i_trial == 0):
      assert approx_equal(es.vectors(), [1,0,0,1])
    assert approx_equal(pai.distance_to_inertia_ellipsoid_surface(
      unit_direction=(1,0)), 6)
    assert pai.distance_to_inertia_ellipsoid_surface(
      unit_direction=(0,0)) == 0
  for i_trial in xrange(10):
    if (i_trial == 0):
      center_of_mass = [0,0]
    else:
      center_of_mass = [rnd(),rnd()]
    points = flex.vec2_double()
    for point in flex.nested_loop([-1,-1], [2,2]):
      points.append((matrix.col([point[0],point[1]*2])
                   + matrix.col(center_of_mass)).elems)
    pai = principal_axes_of_inertia_2d(points=points)
    assert approx_equal(pai.center_of_mass(), center_of_mass)
    assert approx_equal(pai.inertia_tensor(), [24,6,0])
    es = pai.eigensystem()
    assert approx_equal(es.values(), [24,6])
    if (i_trial == 0):
      assert approx_equal(es.vectors(), [1,0,0,1])
  for i_trial in xrange(10):
    if (i_trial == 0):
      center_of_mass = [0,0]
    else:
      center_of_mass = [rnd(),rnd()]
    if (i_trial < 2):
      rot = matrix.sqr([1,0,0,1])
    else:
      rot = euler_angles_as_matrix(
        angles=[random.uniform(0,360),0,0],
        deg=True)
      rot = matrix.sqr([rot.elems[k] for k in [0,1,3,4]])
      theta = random.uniform(0,360)
      csth = math.cos(theta); snth = math.sin(theta)
      rot = matrix.sqr([csth,snth, -snth, csth])
    points = flex.vec2_double()
    for point in [
      [-1,-1],[-1, 1],
      [ 1,-1],[ 1, 1],
      [-1, 0],[ 1, 0],
      [ 0,-1],[ 0, 1],[1,0],[-1,0]]:
      points.append((rot*matrix.col(point)+matrix.col(center_of_mass)).elems)
    pai = principal_axes_of_inertia_2d(points=points)
    assert approx_equal(pai.center_of_mass(), center_of_mass)
    if (i_trial < 2):
      assert approx_equal(pai.inertia_tensor(), [6,8,0])
    es = pai.eigensystem()
    assert approx_equal(es.values(), [8,6])
    if (i_trial < 2):
      assert approx_equal(es.vectors(), [0,1,1,0])
    assert abs(abs(matrix.col(es.vectors()[0:2]).dot(
                   rot*matrix.col([0,1])))-1) < 1.e-3
    assert abs(abs(matrix.col(es.vectors()[2:4]).dot(
                   rot*matrix.col([1,0])))-1) < 1.e-3
    weights = flex.double(points.size(), 1)
    points_plus = points.deep_copy()
    for i_p in [0,2,3,4,5,7,7,7]:
      points_plus.append(points[i_p])
      weights[i_p] += 1
    paip = principal_axes_of_inertia_2d(points=points_plus)
    paiw = principal_axes_of_inertia_2d(points=points, weights=weights)
    assert approx_equal(paip.center_of_mass(), paiw.center_of_mass())
    assert approx_equal(paip.inertia_tensor(), paiw.inertia_tensor())

def explore_inertia_tensor_properties(n_trials=10):
  points = flex.vec3_double([
    (-1,0,0),
    (1,0,0),
    (0,-1,0),
    (0,1,0),
    (0,0,-1),
    (0,0,1)])
  weights = flex.double([2,2,3,3,7,7])
  pai = scitbx.math.principal_axes_of_inertia(
    points=points,
    weights=weights)
  es = pai.eigensystem()
  #
  mt = flex.mersenne_twister(seed=0)
  for i_trial in xrange(n_trials):
    rot_axis = matrix.col(mt.random_double_point_on_sphere())
    rot_angle = 10 + mt.random_double() * 77
    rot_matrix = scitbx.math.r3_rotation_axis_and_angle_as_matrix(
      axis=rot_axis, angle=rot_angle, deg=True)
    #
    rot_points = rot_matrix * points
    rot_pai = scitbx.math.principal_axes_of_inertia(
      points=rot_points,
      weights=weights)
    rot_es = rot_pai.eigensystem()
    #
    c = matrix.sqr(rot_matrix).inverse()
    e = matrix.sym(sym_mat3=rot_pai.inertia_tensor())
    # this proves the transformation law for the inertia tensor
    assert approx_equal(
      (c * e * c.transpose()).as_sym_mat3(),
      pai.inertia_tensor())
    #
    for j_trial in xrange(n_trials):
      v = matrix.col(mt.random_double_point_on_sphere())
      rot_v = matrix.sqr(rot_matrix) * v
      #
      # most intuitive approach (for rwgk)
      def es_distance_to_ellipsoid_surface(es, v):
        assert min(es.values()) > 0
        a,b,c = es.values()
        x,y,z = matrix.sqr(es.vectors()) * v # transform v to eigenvector basis
        # http://mathworld.wolfram.com/Ellipsoid.html
        f = 1/math.sqrt(x*x/(a*a)+y*y/(b*b)+z*z/(c*c))
        return f
      #
      # alternative approach without involving the eigensystem
      def it_distance_to_ellipsoid_surface(inertia_tensor, v):
        iv = matrix.sym(sym_mat3=inertia_tensor).inverse() * v
        return 1/math.sqrt(iv.dot(iv))
      #
      # proves that the intuitive approach works
      d0 = es_distance_to_ellipsoid_surface(es, v)
      d = es_distance_to_ellipsoid_surface(rot_es, rot_v)
      # proves that the alternative approach yields the same results
      assert approx_equal(d, d0)
      d = it_distance_to_ellipsoid_surface(pai.inertia_tensor(), v)
      assert approx_equal(d, d0)
      d = it_distance_to_ellipsoid_surface(rot_pai.inertia_tensor(), rot_v)
      assert approx_equal(d, d0)
      # exercise C++ implementation
      d = pai.distance_to_inertia_ellipsoid_surface(unit_direction=v)
      assert approx_equal(d, d0)
      d = rot_pai.distance_to_inertia_ellipsoid_surface(unit_direction=rot_v)
      assert approx_equal(d, d0)

def exercise_phase_error():
  for deg in [False, True]:
    if (deg): f = 1
    else: f = math.pi/180
    assert approx_equal(signed_phase_error(phi1=-30*f, phi2=270*f, deg=deg),
      -60*f)
    assert approx_equal(signed_phase_error(phi1=330*f, phi2=630*f, deg=deg),
      -60*f)
    assert approx_equal(phase_error(phi1=330*f, phi2=630*f, deg=deg),
      60*f)
    assert approx_equal(nearest_phase(reference=-30*f, other=335*f, deg=deg),
      -25*f)
    assert approx_equal(signed_phase_error(
      phi1=flex.double([-30*f]),
      phi2=flex.double([270*f]), deg=deg), [-60*f])
    assert approx_equal(phase_error(
      phi1=flex.double([-30*f]),
      phi2=flex.double([270*f]), deg=deg), [60*f])
    assert approx_equal(nearest_phase(
      reference=flex.double([-30*f]),
      other=flex.double([345*f]), deg=deg), [-15*f])

def exercise_row_echelon():
  m = flex.int((1,1,1,1))
  m.resize(flex.grid(2,2))
  t = flex.int((2,3))
  t.resize(flex.grid(2,1))
  assert scitbx.math.row_echelon_form_t(m, t) == 1
  assert m.focus() == (1,2)
  assert tuple(m) == (1,1)
  assert tuple(t) == (2,1)
  assert scitbx.math.row_echelon_form(m) == 1
  assert m.focus() == (1,2)
  assert tuple(m) == (1,1)
  m = flex.int((0,-24,0,0,0,-24,24,0,24))
  m.resize(flex.grid(3,3))
  t = flex.int((-3, -6, 0))
  t.resize(flex.grid(3,1))
  assert scitbx.math.row_echelon_form_t(m, t) == 3
  assert tuple(m) == (24,0,24,0,24,0,0,0,24)
  assert tuple(t) == (0,3,6)
  t.resize(flex.grid(3))
  sol = flex.int(3)
  assert scitbx.math.row_echelon_back_substitution_int(m, t, sol) == 8
  assert tuple(sol) == (-2,1,2)
  indep = flex.bool((True,True,True))
  assert scitbx.math.row_echelon_back_substitution_int(
    row_echelon_form=m, independent_flags=indep) == 1
  assert tuple(indep) == (False,False,False)
  #
  for n_cols in xrange(1,5):
    for n_rows in xrange(5):
      for i_trial in xrange(10):
        m = flex.int()
        for i in xrange(n_rows):
          coeffs = flex.int([random.randrange(-5,5) for j in xrange(n_cols)])
          m.extend(coeffs)
        m.resize(flex.grid(n_rows,n_cols))
        rank = scitbx.math.row_echelon_form(m)
        assert m.focus()[0] == rank
        assert m.focus()[1] == n_cols
        indep = flex.bool(n_cols, True)
        scitbx.math.row_echelon_back_substitution_int(
          row_echelon_form=m, independent_flags=indep)
        mm = matrix.rec(m, m.focus())
        s = matrix.col([random.random() for j in xrange(n_cols)])
        sol = flex.double(n_cols, 0)
        sol.set_selected(indep, flex.double(s).select(indep))
        assert scitbx.math.row_echelon_back_substitution_float(
          row_echelon_form=m, solution=sol, v=flex.double(mm * s))
        assert approx_equal(sol, s)
        sol = flex.double(n_cols, 0)
        sol.set_selected(indep, flex.double(s).select(indep))
        assert scitbx.math.row_echelon_back_substitution_float(
          row_echelon_form=m, solution=sol)
        zeros = mm * matrix.col(sol)
        assert approx_equal(zeros, [0]*rank)

def exercise_row_echelon_full_pivoting():
  refp = scitbx.math.row_echelon_full_pivoting
  #
  m = flex.double([
    [ 1,  2,  3,  4,  5,  6],
    [-1, -3,  1,  2, -1,  3],
    [ 2,  1, -1,  3,  4,  2]])
  m_inp = matrix.rec(m, m.all())
  v = [0]*6
  for j in xrange(6):
    for i in xrange(3):
      v[j] += m[i,j]
  e = refp(a_work=m)
  assert list(e.col_perm) == [5, 1, 4, 3, 2, 0]
  assert e.rank == 3
  assert e.nullity == 3
  assert m.all() == (3,6)
  # Is v in the vector space spanned by the rows of m?
  assert e.is_in_row_space(x=flex.double(v), epsilon=1e-15)
  # After that modification, v should not be in that span anymore
  v[2] += 1e-8
  assert not e.is_in_row_space(x=flex.double(v), epsilon=1e-15)
  s = e.back_substitution(free_values=flex.double([1,2,3]))
  assert approx_equal(s, [
    3.0, -0.42857142857142849, 2.0,
    1.0, -1.0816326530612246, -1.1224489795918366])
  assert approx_equal(m_inp * matrix.col(s), [0,0,0])
  #
  # Let's test with a row rank deficient matrix m now
  m = flex.double([
    [ 1,  2,  3,  4,  5,  6],
    [-1, -3,  1,  2, -1,  3],
    [ 2,  1,  7, 10,  9, 15]])
  m_inp = matrix.rec(m, m.all())
  v = [0]*6
  for j in xrange(6):
    m[2,j] =   m[1,j] + 2*m[0,j]
    v[j]   = 2*m[1,j] +   m[0,j]
  e = refp(a_work=m, min_abs_pivot=1e-15)
  assert list(e.col_perm) == [5, 1, 2, 3, 4, 0]
  assert e.rank == 2
  assert e.nullity == 4
  assert e.is_in_row_space(x=flex.double(v), epsilon=1e-15)
  v[4] += 1e-9
  assert not e.is_in_row_space(x=flex.double(v), epsilon=1e-15)
  s = e.back_substitution(free_values=flex.double([-3,1,4,-2]))
  assert approx_equal(s, [-2.0, -2.3749999999999991, -3.0, 1.0, 4.0, -1.375])
  assert approx_equal(m_inp * matrix.col(s), [0,0,-2])
  #
  try: refp(a_work=flex.double())
  except RuntimeError, e:
    assert str(e) == "a_work matrix must be two-dimensional."
  else: raise Exception_expected
  for v in [0,1,-1]:
    for nr in xrange(5):
      for nc in xrange(5):
        a = flex.double(flex.grid(nr, nc), v)
        e = refp(a_work=a)
        assert e.rank == min(abs(v), nr, nc)
        assert e.nullity == nc - e.rank
        a = flex.double(flex.grid(nr, nc), v)
        b = flex.double(nr)
        e = refp(a_work=a, b_work=b)
        assert e.rank == min(abs(v), nr, nc)
        assert e.nullity == nc - e.rank
  #
  mt = flex.mersenne_twister(seed=0)
  for i_trial in xrange(10):
    for nr in xrange(1,5):
      for nc in xrange(1,5):
        a = mt.random_double(size=nr*nc)*2-1
        a.reshape(flex.grid(nr, nc))
        x = mt.random_double(size=nc)*2-1
        b = a.matrix_multiply(x)
        aw = a.deep_copy()
        bw = b.deep_copy()
        e = refp(a_work=aw, b_work=bw)
        assert e.rank == min(nr, nc) # assumes (random) linear indepence
        assert e.nullity == nc - e.rank
        for v in [0,1,-1]:
          ex = e.back_substitution(
            free_values=flex.double(e.nullity, v),
            epsilon=1e-10)
          assert ex is not None
          assert ex.size() == nc
          eb = a.matrix_multiply(ex)
          assert approx_equal(eb, b)
        ex = e.back_substitution(
          free_values=x.select(flex.size_t(iter(e.col_perm[e.rank:]))),
          epsilon=1e-10)
        assert approx_equal(ex, x)
  #
  a = flex.double([
    [1, 2,  ],
    [1, 1.99],
    [1, 1.98 ]])
  aw = a.deep_copy()
  e = refp(a_work=aw)
  assert approx_equal(aw, [2, 1, 0, 0.01, 0, 0])
  assert e.rank == 2
  aw = a.deep_copy()
  e = refp(a_work=aw, min_abs_pivot=0.1)
  assert e.rank == 1
  assert approx_equal(aw, [2, 1, 0, 0.005, 0, 0.01])
  aw_rank_1 = aw
  aw = a.deep_copy()
  e = refp(a_work=aw, max_rank=1)
  assert e.rank == 1
  assert approx_equal(aw, aw_rank_1)
  aw = a.deep_copy()
  e = refp(a_work=aw, max_rank=0)
  assert e.rank == 0
  assert approx_equal(aw, a)
  #
  n_no_solution_with_epsilon_zero = 0
  for singular_a,rank in [
        ((3,0,0,0,-2,0,0,0,0), 2),
        ((0,0,3,0,0,0,0,0,0), 1),
        ((0,0,0,1e-15,0,0,0,0,0), 0)]:
    for i_trial in xrange(10):
      r = matrix.sqr(mt.random_double_r3_rotation_matrix())
      a = flex.double(r * matrix.sqr(singular_a) * r.transpose())
      a.reshape(flex.grid(3,3))
      assert approx_equal(a.matrix_determinant_via_lu(), 0, eps=1e-12)
      x = flex.double([4.3,-2.1,8.2])
      b = a.matrix_multiply(x)
      aw = a.deep_copy()
      bw = b.deep_copy()
      e = refp(a_work=aw, b_work=bw, min_abs_pivot=1e-12)
      assert e.rank == rank
      ex = e.back_substitution(
        free_values=flex.double(e.nullity, 9.3),
        epsilon=1e-12)
      assert ex is not None
      eb = a.matrix_multiply(ex)
      assert approx_equal(eb, b)
      ex = e.back_substitution(free_values=flex.double(e.nullity, 0))
      if (ex is None): n_no_solution_with_epsilon_zero += 1
  assert n_no_solution_with_epsilon_zero > 20
  #
  # http://www.mathworks.com/access/helpdesk/help/techdoc/ref/pinv.html
  #   2008-12-01
  a = flex.double([float(v) for v in """
    64     2     3    61    60     6
     9    55    54    12    13    51
    17    47    46    20    21    43
    40    26    27    37    36    30
    32    34    35    29    28    38
    41    23    22    44    45    19
    49    15    14    52    53    11
     8    58    59     5     4    62""".split()])
  a.reshape(flex.grid(8,6))
  b = flex.double(8, 260)
  aw = a.deep_copy()
  bw = b.deep_copy()
  e = scitbx.math.row_echelon_full_pivoting(
    a_work=aw, b_work=bw, min_abs_pivot=1e-12)
  assert e.rank == 3
  x = e.back_substitution(free_values=flex.double(e.nullity), epsilon=1e-12)
  assert approx_equal(a.matrix_multiply(x), b)
  assert approx_equal(sorted(x), [-1,0,0,0,4,5])

def exercise_solve_a_x_eq_b_min_norm_given_a_sym_b_col():
  def girs(a, relative_min_abs_pivot=1e-12, absolute_min_abs_pivot=0):
    es = scitbx.linalg.eigensystem.real_symmetric(
      m=a,
      relative_epsilon=relative_min_abs_pivot,
      absolute_epsilon=absolute_min_abs_pivot)
    return es.generalized_inverse_as_packed_u().matrix_packed_u_as_symmetric()
  mt = flex.mersenne_twister(seed=0)
  for bits in xrange(8):
    d = [1.23, 2.34, 0.58]
    x = [-0.19, -0.44, 0.83]
    if (bits    % 2): d[0] = x[0] = 0
    if (bits//2 % 2): d[1] = x[1] = 0
    if (bits//4 % 2): d[2] = x[2] = 0
    a = matrix.diag(d)
    x = matrix.col(x)
    b = a * x
    xs = scitbx.math.solve_a_x_eq_b_min_norm_given_a_sym_b_col(a=a, b=b)
    assert approx_equal(xs, x)
    for i_trial in xrange(10):
      if (i_trial == 0):
        r = matrix.identity(n=3)
      else:
        r = matrix.sqr(mt.random_double_r3_rotation_matrix())
      ar = r * a * r.transpose()
      xr = r * x
      br = r * b
      assert approx_equal(ar * xr, br)
      xs = scitbx.math.solve_a_x_eq_b_min_norm_given_a_sym_b_col(a=ar, b=br)
      assert approx_equal(xs, xr)
      #
      ari = matrix.sqr(girs(a=ar.as_flex_double_matrix()))
      assert approx_equal(ari * br, xr)
      #
      if (tntbx is not None):
        arit = tntbx.generalized_inverse(ar.as_flex_double_matrix())
        xs = matrix.sqr(arit) * br
        assert approx_equal(xs, xr)
        assert approx_equal(ari, arit)
  #
  a = flex.double([[1e-15]])
  b = flex.double([1e-14])
  x = scitbx.math.solve_a_x_eq_b_min_norm_given_a_sym_b_col(a=a, b=b)
  assert approx_equal(x, [10])
  ai = matrix.sqr(girs(a=a))
  x = ai * matrix.col(b)
  assert approx_equal(x, [10])
  assert a[0] == 1e-15
  assert b[0] == 1e-14
  x = scitbx.math.solve_a_x_eq_b_min_norm_given_a_sym_b_col(a=a, b=b,
    absolute_min_abs_pivot=1e-12)
  assert x[0] == 0
  ai = matrix.sqr(girs(a=a, absolute_min_abs_pivot=1e-12))
  x = ai * matrix.col(b)
  assert x[0] == 0
  b[0] = 1e-10
  x = scitbx.math.solve_a_x_eq_b_min_norm_given_a_sym_b_col(a=a, b=b,
    absolute_min_abs_pivot=1e-12)
  assert x is None
  #
  def compare(a, n_trials):
    for i_trial in xrange(n_trials):
      if (i_trial == 0):
        ar = a
      else:
        r = matrix.sqr(mt.random_double_r3_rotation_matrix())
        ar = r * a * r.transpose()
      ari = matrix.sqr(girs(
        a=ar.as_flex_double_matrix(), absolute_min_abs_pivot=1.e-12))
      if (tntbx is not None):
        arit = matrix.sqr(
          tntbx.generalized_inverse(ar.as_flex_double_matrix()))
        mismatch = (ari-arit).norm_sq() / max(1, max([abs(e) for e in ari]))
        if (mismatch > 1e-10):
          print ar.elems
          print ari.elems
          print arit.elems
          raise AssertionError, mismatch
  for i_trial in xrange(10):
    x,y,z = flex.random_double(size=3)*2-1
    a = matrix.sqr([
      x,y,x,
      y,z,y,
      x,y,x])
    compare(a, n_trials=10)

def exercise_tensor_rank_2():
  g = (2,3,5,0.2,0.3,0.5)
  assert approx_equal(scitbx.math.tensor_rank_2_gradient_transform(
    a=(1,0,0,0,1,0,0,0,1), g=g), g)
  a = (-0.00266542,0.386546, 0.22833,
        0.263694, -0.660647, 0.896465,
        0.888726, -0.996946,-0.521507)
  assert approx_equal(matrix.sqr(a).determinant(), 0.431857368657)
  ga = scitbx.math.tensor_rank_2_gradient_transform(a=a, g=g)
  assert approx_equal(ga,
    [4.2741119386687805, 6.7403365850628001, 3.6465242980395001,
     -10.209907479357136, -2.8163934767020788, 1.6344744599549008])
  gtmx = scitbx.math.tensor_rank_2_gradient_transform_matrix(a=a)
  assert gtmx.focus() == (6,6)
  assert approx_equal(gtmx.matrix_multiply(flex.double(g)), ga)

def exercise_minimum_covering_sphere(epsilon=1.e-3):
  s3 = sphere_3d(center=[1,2,3], radius=4)
  assert approx_equal(s3.center(), [1,2,3])
  assert approx_equal(s3.radius(), 4)
  s3 = s3.expand(additional_radius=2)
  assert approx_equal(s3.center(), [1,2,3])
  assert approx_equal(s3.radius(), 6)
  assert approx_equal(
    s3.expand_relative(additional_relative_radius=0.1).radius(), 6.6)
  assert s3.is_inside(point=[1,2,3])
  assert s3.is_inside(point=[1,2,3+6-1.e-6])
  assert not s3.is_inside([1,2,3+6+1.e-6])
  assert approx_equal(s3.box_min(), [1-6,2-6,3-6])
  assert approx_equal(s3.box_max(), [1+6,2+6,3+6])
  for i_impl,mcs_impl in enumerate([scitbx.math.minimum_covering_sphere_3d,
                                    scitbx.math.minimum_covering_sphere_nd]):
    def wrap_points(points):
      if (i_impl == 0): return flex.vec3_double(points)
      return [matrix.col(point) for point in points]
    points = wrap_points([])
    mcs = mcs_impl(points=points)
    assert mcs.n_iterations() == 0
    assert approx_equal(mcs.center(), (0,0,0))
    assert approx_equal(mcs.radius(), 1)
    if (i_impl == 0):
      assert mcs.is_inside(mcs.center()) # base class method
    mcs = mcs_impl(
      points=points,
      epsilon=1.e-6,
      radius_if_one_or_no_points=3,
      center_if_no_points=(2,3,5))
    assert mcs.n_iterations() == 0
    assert approx_equal(mcs.center(), (2,3,5))
    assert approx_equal(mcs.radius(), 3)
    points = wrap_points([(3,4,5)])
    mcs = mcs_impl(
      points=points,
      epsilon=1.e-6,
      radius_if_one_or_no_points=5,
      center_if_no_points=(2,3,5))
    assert mcs.n_iterations() == 0
    assert approx_equal(mcs.center(), (3,4,5))
    assert approx_equal(mcs.radius(), 5)
    points = wrap_points([(0,0,0),(1,0,0),(0,1,0),(1,1,1)])
    mcs = mcs_impl(points=points)
    assert mcs.n_iterations() > 0
    assert approx_equal(mcs.center(), (0.5,0.5,0.5), eps=1.e-3)
    assert approx_equal(mcs.radius(), math.sqrt(3)/2, eps=1.e-5)
  eps = epsilon*10
  eps_loose = eps*10
  for i,j,k in flex.nested_loop((1,1,1),(2,3,2),False):
    for shift in [(0,0,0),(2,3,4),(-3,-5,2)]:
      for poly_index in xrange(1,2):
        if (poly_index == 0):
          # cube
          points = flex.vec3_double(
            [(matrix.col(t)+matrix.col(shift)).elems for t in [
            (0,0,0),
            (0,0,k),
            (0,j,0),
            (0,j,k),
            (i,0,0),
            (i,0,k),
            (i,j,0),
            (i,j,k)]])
          expected_center = (matrix.col(shift) + matrix.col([i,j,k])/2.).elems
          expected_radius = math.sqrt(i**2+j**2+k**2)/2
        else:
          # tetrahedron
          z = 1/math.sqrt(2)*k
          points = flex.vec3_double(
            [(matrix.col(t)/2.+matrix.col(shift)).elems for t in [
            (-i,0,z),
            (i,0,z),
            (0,-j,-z),
            (0,j,-z)]])
          if (i == j and j == k):
            expected_center = shift
            expected_radius = max(
              [abs(matrix.col(points[0])-matrix.col(shift))
                for point in points])
          else:
            expected_center = None
            expected_radius = None
        mcs = minimum_covering_sphere(points, epsilon=epsilon)
        if (expected_center is None):
          expected_center = mcs.center()
          expected_radius = mcs.radius()
        assert approx_equal(mcs.center(), expected_center, eps=eps)
        assert approx_equal(mcs.radius(), expected_radius, eps=eps)
        if (poly_index == 0):
          assert mcs.n_iterations() == 0
        points.append(expected_center)
        mcs = minimum_covering_sphere(points, epsilon=epsilon)
        assert approx_equal(mcs.center(), expected_center, eps=eps)
        assert approx_equal(mcs.radius(), expected_radius, eps=eps)
        if (poly_index == 0):
          assert mcs.n_iterations() <= 1
        r = random.random
        for i_addl in xrange(3):
          points.append(
            (matrix.col(expected_center)
             + matrix.col([r(),r(),r()]).normalize()*expected_radius).elems)
          mcs = minimum_covering_sphere(points, epsilon=epsilon)
          assert approx_equal(mcs.center(), expected_center, eps=eps_loose)
          assert approx_equal(mcs.radius(), expected_radius, eps=eps)
        # also exercise the Python implementation
        mcs = minimum_covering_sphere(
          points=[matrix.col(point) for point in points],
          epsilon=epsilon)
        assert approx_equal(mcs.center(), expected_center, eps=eps_loose)
        assert approx_equal(mcs.radius(), expected_radius, eps=eps)
  # exercise Python implementation with sets of 2-dimensional points
  for i,j in flex.nested_loop((1,1),(2,3),False):
    for shift in [(0,0),(3,4),(-3,-5)]:
      # square
      points = [matrix.col(t)+matrix.col(shift) for t in [
        (0,0),
        (0,j),
        (i,0),
        (i,j)]]
      expected_center = (matrix.col(shift) + matrix.col([i,j])/2.).elems
      expected_radius = math.sqrt(i**2+j**2)/2
      mcs = minimum_covering_sphere(points, epsilon=epsilon)
      if (expected_center is None):
        expected_center = mcs.center()
        expected_radius = mcs.radius()
      assert approx_equal(mcs.center(), expected_center, eps=eps)
      assert approx_equal(mcs.radius(), expected_radius, eps=eps)
      assert mcs.n_iterations() == 0
      points.append(matrix.col(expected_center))
      mcs = minimum_covering_sphere(points, epsilon=epsilon)
      assert approx_equal(mcs.center(), expected_center, eps=eps)
      assert approx_equal(mcs.radius(), expected_radius, eps=eps)
      assert mcs.n_iterations() <= 1
      r = random.random
      for i_addl in xrange(3):
        points.append(
          matrix.col(expected_center)
          + matrix.col([r(),r()]).normalize()*expected_radius)
        mcs = minimum_covering_sphere(points, epsilon=epsilon)
        assert approx_equal(mcs.center(), expected_center, eps=eps_loose)
        assert approx_equal(mcs.radius(), expected_radius, eps=eps)

def exercise_icosahedron():
  ico = icosahedron(level=0)
  for level in xrange(6):
    ico = icosahedron(level=level)
    assert ico.level == level
    if (level == 0):
      assert ico.sites.size() == 12
    else:
      assert ico.sites.size() == 80 * 4**(level-1)
    assert approx_equal(ico.sites.mean(), [0,0,0])
    assert approx_equal(ico.sites.dot(), [1]*ico.sites.size())
    d = ico.next_neighbors_distance()
    m = flex.min((ico.sites[1:] - ico.sites[0]).dot())**0.5
    if (level == 0):
      assert approx_equal(d, m)
    else:
      assert d > m
      assert d/2 < m

def exercise_basic_statistics():
  x = flex.double([])
  s = scitbx.math.basic_statistics(values=x)
  assert s.n == 0
  assert approx_equal(s.min, -1)
  assert approx_equal(s.max, -1)
  assert approx_equal(s.max_absolute, -1)
  assert approx_equal(s.sum, -1)
  assert approx_equal(s.mean, -1)
  assert approx_equal(s.mean_absolute_deviation_from_mean, -1)
  assert approx_equal(s.biased_variance, -1)
  assert approx_equal(s.biased_standard_deviation, -1)
  assert approx_equal(s.bias_corrected_variance, -1)
  assert approx_equal(s.bias_corrected_standard_deviation, -1)
  assert approx_equal(s.skew, -1)
  assert approx_equal(s.kurtosis, -1)
  assert approx_equal(s.kurtosis_excess, -1)
  x = flex.double([-7])
  s = scitbx.math.basic_statistics(values=x)
  assert s.n == 1
  assert approx_equal(s.min, -7)
  assert approx_equal(s.max, -7)
  assert approx_equal(s.max_absolute, 7)
  assert approx_equal(s.sum, -7)
  assert approx_equal(s.mean, -7)
  assert approx_equal(s.mean_absolute_deviation_from_mean, 0)
  assert approx_equal(s.biased_variance, 0)
  assert approx_equal(s.biased_standard_deviation, 0)
  assert approx_equal(s.bias_corrected_variance, -1)
  assert approx_equal(s.bias_corrected_standard_deviation, -1)
  assert approx_equal(s.skew, -1)
  assert approx_equal(s.kurtosis, -1)
  assert approx_equal(s.kurtosis_excess, -1)
  x = flex.double([1,2,3,4,5])
  s = scitbx.math.basic_statistics(values=x)
  assert s.n == 5
  assert approx_equal(s.min, 1)
  assert approx_equal(s.max, 5)
  assert approx_equal(s.max_absolute, 5)
  assert approx_equal(s.sum, 15)
  assert approx_equal(s.mean, 3)
  assert approx_equal(s.mean_absolute_deviation_from_mean, 1.2)
  assert approx_equal(s.biased_variance, 2)
  assert approx_equal(s.biased_standard_deviation, math.sqrt(2))
  assert approx_equal(s.bias_corrected_variance, 2.5)
  assert approx_equal(s.bias_corrected_standard_deviation, math.sqrt(2.5))
  assert approx_equal(s.skew, 0)
  assert approx_equal(s.kurtosis, 1.7)
  assert approx_equal(s.kurtosis_excess, -1.3)
  x = flex.double([1,1,1])
  s = scitbx.math.basic_statistics(values=x)
  assert s.n == 3
  assert approx_equal(s.min, 1)
  assert approx_equal(s.max, 1)
  assert approx_equal(s.max_absolute, 1)
  assert approx_equal(s.sum, 3)
  assert approx_equal(s.mean, 1)
  assert approx_equal(s.mean_absolute_deviation_from_mean, 0)
  assert approx_equal(s.biased_variance, 0)
  assert approx_equal(s.biased_standard_deviation, 0)
  assert approx_equal(s.bias_corrected_variance, 0)
  assert approx_equal(s.bias_corrected_standard_deviation, math.sqrt(0))
  assert approx_equal(s.skew, -1)
  assert approx_equal(s.kurtosis, -1)
  assert approx_equal(s.kurtosis_excess, -1)
  f = StringIO()
  s.show(f=f)
  assert len(f.getvalue().splitlines()) == 14
  for i_trial in xrange(10):
    x = flex.random_double(size=2+int(random.random()*10))
    s = scitbx.math.basic_statistics(values=x)
    assert s.n == x.size()
    assert approx_equal(s.min, flex.min(x))
    assert approx_equal(s.max, flex.max(x))
    assert approx_equal(s.max_absolute, max(-flex.min(x), flex.max(x)))
    assert approx_equal(s.sum, flex.sum(x))
    assert approx_equal(s.mean, flex.mean(x))
    d = x-flex.mean(x)
    assert approx_equal(s.mean_absolute_deviation_from_mean,
      flex.mean(flex.abs(d)))
    assert approx_equal(s.biased_variance, flex.sum(d*d) / s.n)
    assert approx_equal(s.biased_standard_deviation,
      math.sqrt(s.biased_variance))
    assert approx_equal(s.bias_corrected_variance,
      flex.sum(flex.pow2(d)) / (s.n-1))
    assert approx_equal(s.bias_corrected_standard_deviation,
      math.sqrt(s.bias_corrected_variance))
    assert approx_equal(s.skew,
      (flex.sum(d*d*d)/s.n) / (flex.sum(d*d)/s.n)**(3/2.))
    assert approx_equal(s.kurtosis,
      (flex.sum(d*d*d*d)/s.n) / (flex.sum(d*d)/s.n)**2)
    assert approx_equal(s.kurtosis_excess, s.kurtosis-3)

def exercise_median():
  from scitbx.math import median_statistics

  try:
    median_statistics(flex.double())
  except RuntimeError:
    pass
  else:
    raise Exception_expected

  stats = median_statistics(flex.double((1,)))
  assert stats.median == 1
  assert stats.median_absolute_deviation == 0

  for i in xrange(5):
    stats = median_statistics(flex.double((5, 1)))
    assert stats.median == 3
    assert stats.median_absolute_deviation == 2

  for i in xrange(5):
    stats = median_statistics(flex.double((5, 1, 2)))
    assert stats.median == 2
    assert stats.median_absolute_deviation == 1

  data = flex.double((1, 1, 2, 2, 4, 6, 9))
  for i in xrange(10):
    data_ = data.select(flex.random_permutation(len(data)))
    stats = median_statistics(data_)
    assert stats.median == 2
    assert stats.median_absolute_deviation == 1

  data = flex.double((1, 1, 2, 4, 6, 9))
  for i in xrange(10):
    data_ = data.select(flex.random_permutation(len(data)))
    stats = median_statistics(data_)
    assert stats.median == 3
    assert stats.median_absolute_deviation == 2

def exercise_slatec_dlngam():
  def cmp(a, b):
    if (abs(a) < 1):
      assert approx_equal(a, b, eps=1.e-10)
    else:
      assert approx_equal((a-b)/(abs(a+b)), 0, eps=1.e-10)
  try: slatec_dlngam(x=0)
  except RuntimeError, e:
    assert str(e)=="slatec: dgamma: x is 0 (nerr=4, level=2)"
  else: raise Exception_expected
  try: slatec_dlngam(x=-1)
  except RuntimeError, e:
    assert str(e)=="slatec: dgamma: x is a negative integer (nerr=4, level=2)"
  else: raise Exception_expected
  for i in xrange(1,10000):
    x = i/100.
    cmp(slatec_dgamma(x=x), gamma_complete(x))
  try: slatec_dlngam(x=0)
  except RuntimeError, e:
    assert str(e)=="slatec: dgamma: x is 0 (nerr=4, level=2)"
  else: raise Exception_expected
  try: slatec_dlngam(-1)
  except RuntimeError, e:
    assert str(e)=="slatec: dgamma: x is a negative integer (nerr=4, level=2)"
  else: raise Exception_expected
  assert approx_equal(slatec_dlngam(-1+1.e-8), 18.4206807543)
  assert approx_equal(slatec_dlngam(-1-1.e-8), 18.4206807458)
  assert eps_eq(slatec_dlngam( 2.53273727e+305),  1.77853307723e+308)
  try: slatec_dlngam(-2.53273727e+305)
  except RuntimeError, e:
    assert str(e)=="slatec: dlngam: x is a negative integer (nerr=3, level=2)"
  else: raise Exception_expected
  for x in [2.53273728e+305, -2.53273728e+305]:
    try: slatec_dlngam(x=x)
    except RuntimeError, e:
      assert str(e) == \
        "slatec: dlngam: abs(x) so big dlngam overflows (nerr=2, level=2)"
    else: raise Exception_expected
  for x,y in [
        (-0.9, 2.35807316739203), (-0.8, 1.74720737374499),
        (-0.7, 1.45247293875681), (-0.6, 1.30750344146777),
        (-0.5, 1.26551212348465), (-0.4, 1.31452458994339),
        (-0.3, 1.4648400508576), (-0.2, 1.76149759083394),
        (-0.1, 2.36896133272879), (0.1, 2.25271265173421),
        (0.2, 1.52406382243078), (0.3, 1.09579799481808),
        (0.4, 0.796677817701784), (0.5, 0.5723649429247),
        (0.6, 0.398233858069235), (0.7, 0.260867246531667),
        (0.8, 0.152059678399838), (0.9, 0.066376239734743),
        (-95.7, -342.652344377166), (-95.4, -341.444636043021),
        (-94.4, -336.886557464567), (-94.3, -336.269573770848),
        (-89.7, -315.444865680315), (-85.3, -295.745018427395),
        (-81.4, -278.633885462703), (-75.9, -253.468382259545),
        (-70.3, -230.359702035999), (-61.4, -193.193052861824),
        (-61.2, -191.887058269308), (-56.8, -173.90971361753),
        (-54.9, -165.606807044062), (-53.4, -160.729586596082),
        (-52.1, -154.437925415), (-45.4, -129.457867985401),
        (-42.5, -118.504844660495), (-28.9, -68.5996759295182),
        (-28.5, -68.4243510349742), (8.7, 9.96776168512864),
        (14.4, 23.5991967127359), (48.1, 137.188902640497),
        (52.4, 153.98778093456), (57.7, 175.181093095627),
        (58.7, 179.236350269141), (76.1, 252.322882401268),
        (80.1, 269.728736878324), (91.7, 321.309088278786),
        (95.2, 337.171114368332), (97.3, 346.750737141662),
        (99.2, 355.457300594627),
        (-1.25992104989487, 1.33050808569476),
        (3.1748021039364, 0.860352839192692),
        (-12.6992084157456, -20.4177393801341),
        (-20.158736798318, -40.9333927025327),
        (-322.539788773088, -1543.17783526817),
        (-812.749338607718, -4635.79483859481),
        (134217728, 2377663536.65922),
        (-213057362.619982, -3872759028.80297),
        (-426114725.239963, -8040878268.75138),
        (-1352829926.21012, -27091047640.5208),
        (-1704458900.95985, -34526394809.1753),
        (6817835603.83941, 147557106245.627),
        (-10822639409.6809, -239233427036.051),
        (-1385297844439.16, -37343385724961),
        (-3490731829165.78, -97325556732394.7),
        (-22164765511026.5, -658947950486145),
        (-354636248176425, -1.15264276699695e+16),
        (2.90518014506127e+18, 1.20602822017773e+20),
        (4.64828823209803e+19, 2.05852306758472e+21),
        (1.18059162071741e+21, 5.61020711097905e+22),
        (5.94980893708548e+21, 2.92359605678122e+23),
        (1.54742504910673e+26, 9.17681929135993e+27),
        (6.43926366825732e+31, 4.65188840905571e+33),
        (1.72852667909289e+40, 1.58420633670403e+42),
        (3.95432630437266e+57, 5.20476971689902e+59),
        (3.16346104349813e+58, 4.22959809661257e+60),
        (6.37713785476111e+59, 8.71788160127098e+61),
        (1.27542757095222e+60, 1.75241691050233e+62),
        (4.0492301356776e+60, 5.61035543530293e+62),
        (3.36999333339383e+66, 5.12864211104658e+68),
        (2.69599466671506e+67, 4.15897532189841e+69),
        (5.47791970565651e+69, 8.74161479134395e+71),
        (6.90174634679056e+69, 1.10296909057573e+72),
        (2.76069853871623e+70, 4.45014777047268e+72),
        (2.80469488929613e+72, 4.65067274510429e+74),
        (7.29444453207001e+76, 1.28370246361623e+79),
        (1.1671111251312e+78, 2.08628313321737e+80),
        (1.18571099379012e+80, 2.1743230279922e+82),
        (2.08592483976651e+93, 4.46128480838909e+95),
        (5.42506890849806e+97, 1.21544211995612e+100),
        (1.92392608380832e+112, 4.95495423759051e+114),
        (3.10271356343285e+114, 8.14856285749288e+116),
        (1.24108542537314e+115, 3.27663024026578e+117),
        (3.15216049571156e+116, 8.42408458530677e+118),
        (3.20239544759368e+118, 8.70631856055695e+120),
        (2.56191635807495e+119, 7.01832840145846e+121),
        (1.63962646916797e+121, 4.55992032478993e+123),
        (3.33151332094993e+123, 9.44222146520736e+125),
        (4.33229639706377e+127, 1.26890537358424e+130),
        (1.09167028500887e+128, 3.20753035674239e+130),
        (2.53719906956893e+147, 8.58616250449499e+149),
        (4.05951851131028e+148, 1.38504137596138e+151),
        (1.02293456496754e+149, 3.49953949189293e+151),
        (2.61871248631691e+151, 9.10403335337293e+153),
        (6.65111512893076e+152, 2.3337946141639e+155),
        (3.54267261962962e+160, 1.3061072142834e+163),
        (4.02035240429694e+176, 1.63084126214792e+179),
        (2.55276465434554e+177, 1.04023816504482e+180),
        (1.39234637988959e+188, 6.01795758359681e+190),
        (1.13162859936811e+191, 4.96691494590327e+193),
        (3.62121151797796e+192, 1.60196294545853e+195),
        (3.70812059440943e+195, 1.66611278950145e+198),
        (1.97510598530428e+203, 9.22582663285365e+205),
        (5.51565226310199e+216, 2.74715869462202e+219),
        (2.24142136441315e+219, 1.12984054184522e+222),
        (1.85074578797902e+224, 9.53864476740804e+226),
        (4.81341572835509e+228, 2.52974266023021e+231),
        (7.94889263257963e+233, 4.27312881295699e+236),
        (2.56383415069219e+236, 1.39306335150614e+239),
        (1.35485608003746e+243, 7.57136760625491e+245),
        (2.97936002792839e+255, 1.74963187938126e+258),
        (3.09948530198153e+260, 1.8559822138551e+263),
        (3.99882933842564e+263, 2.42315414346892e+266),
        (3.4622310392507e+274, 2.18518791213377e+277),
        (2.97403381695557e+284, 1.94508974776453e+287),
        (3.80676328570312e+286, 2.50818540780841e+289),
        (9.9792015476736e+291, 6.69956455295144e+294),
        (6.96694329442493e+304, 4.88331897497815e+307)]:
    cmp(y, slatec_dlngam(x=x))
  cmath_lgamma = getattr(scitbx.math, "cmath_lgamma", None)
  if (cmath_lgamma is not None):
    print "Testing compatibility of cmath_lgamma and slatec_dlngam...",
    for i in xrange(-1000,1000):
      if (i <= 0 and i % 10 == 0): continue
      x = i/10.
      assert approx_equal(slatec_dlngam(x), cmath_lgamma(x), eps=1.e-10)
    cmath_lgamma_max_x = 5.e15 # larger values lead to floating-point
                               # exceptions on some platforms
    v = 2**(1/3.)
    x = v
    while True:
      try: s = slatec_dlngam(x)
      except RuntimeError, e:
        assert str(e) == \
          "slatec: dlngam: abs(x) so big dlngam overflows (nerr=2, level=2)"
        break
      if (x < cmath_lgamma_max_x):
        m = cmath_lgamma(x)
        cmp(s, m)
      try: s = slatec_dlngam(-x)
      except RuntimeError, e:
        assert str(e) in [
          "slatec: dlngam: x is a negative integer (nerr=3, level=2)",
          "slatec: dgamma: x is a negative integer (nerr=4, level=2)"]
      else:
        if (x < cmath_lgamma_max_x):
          m = cmath_lgamma(-x)
          cmp(s, m)
      x *= v
    print "OK"

def exercise_slatec_dbinom():
  f = scitbx.math.slatec_dlnrel
  try: f(-1)
  except RuntimeError, e:
    assert str(e) == \
      "slatec: dlnrel: x is le -1 (nerr=2, level=2)"
  else: raise Exception_expected
  assert approx_equal(f(-1+1.e-10), -23.0258508472)
  assert approx_equal(f(0.374), 0.3177261938)
  assert approx_equal(f(0.376), 0.319180739511)
  assert eps_eq(f(-0.4), -0.510825623766)
  assert eps_eq(f(0), 0.0)
  assert eps_eq(f(0.3), 0.262364264467)
  assert eps_eq(f(0.4), 0.336472236621)
  f = scitbx.math.slatec_dbinom
  try: f(n=0, m=1)
  except RuntimeError, e:
    assert str(e) == "slatec: dbinom: n lt m (nerr=2, level=2)"
  else: raise Exception_expected
  expected = [
    1, 2, 1, 3, 3, 1, 4, 6, 4, 1, 5, 10, 10, 5, 1, 6, 15, 20, 15, 6,
    1, 7, 21, 35, 35, 21, 7, 1, 8, 28, 56, 70, 56, 28, 8, 1, 9, 36, 84,
    126, 126, 84, 36, 9, 1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1]
  i = 0
  for n in xrange(1,11):
    for m in xrange(1,n+1):
      assert approx_equal(f(n=n, m=m), expected[i])
      i += 1
  assert eps_eq(f(100, 10), 1.73103095E+13)
  assert eps_eq(f(100, 33), 2.94692427E+26)
  assert eps_eq(f(1000, 100), 6.38505119E+139)
  assert eps_eq(f(1000, 333), 5.77613455E+274)
  nms = [
    (5, 2), (9, 6), (8, 3), (9, 1), (8, 2),
    (6, 1), (8, 4), (7, 5), (7, 3), (8, 7),
    (93, 70), (64, 57), (76, 66), (55, 22), (70, 2),
    (90, 85), (78, 4), (82, 19), (99, 6), (71, 5),
    (957, 516), (896, 665), (909, 253), (579, 74), (653, 651),
    (820, 581), (638, 290), (697, 533), (937, 695), (725, 78)]
  expected = [
    10, 84, 56, 9, 28, 6, 70, 21, 35, 8,
    3.73549788E+21, 621216192., 9.54526729E+11, 1.30085363E+15, 2415.,
    43949268., 1426425., 1.97103824E+18, 1.12052926E+09, 13019909.,
    1.66252414E+285, 3.80970836E+220, 8.43685887E+231, 6.3253529E+94, 212878.,
    2.45633786E+213, 2.58081251E+189, 5.06707246E+163, 8.53013061E+230,
    1.53361301E+106]
  for nm,e in zip(nms,expected):
    assert eps_eq(f(*nm), e)
  assert eps_eq(f(n=2**32-1,m=2**5), 6.83193552992e+272)
  try: f(n=2**32-1,m=2**6)
  except RuntimeError, e:
    assert str(e) == \
      "slatec: dbinom: result overflows" \
      " because n and/or m too big (nerr=3, level=2)"
  else: raise Exception_expected

def exercise_unimodular_generator(forever):
  ug = scitbx.math.unimodular_generator
  g = ug(range=0)
  assert g.at_end()
  g = ug(range=1)
  assert not g.at_end()
  n = 0
  while (not g.at_end()):
    assert matrix.rec(g.next(), (3,3)).determinant() == 1
    n += 1
  assert n == 3480
  assert ug(range=0).count() == 0
  assert ug(range=1).count() == 3480
  assert ug(range=2).count() == 67704
  assert ug(range=3).count() == 640824
  assert len(list(ug(range=1).all())) == 3480
  for range in count():
    timer = user_plus_sys_time()
    n = ug(range=range).count()
    print "unimodular range %d: count=%d, time=%.2f s" % (
      range, n, timer.elapsed())
    if (range == 4 and not forever):
      break

def exercise_least_squares_plane():
  points = [ matrix.col(x) for x in [(1, 2, 3), (-1, -2, -3), (1, 1, 1),
                                     (1.2, 2.1, 2.9), (-0.9, -2.1, -3.1),
                                     (1.1, 0.8, 1.2)] ]
  def distance(n,d):
    u = n/d
    return sum([ (u.dot(x) - 1)**2 for x in points ])
  flex_points = flex.vec3_double(points)
  p = scitbx.math.least_squares_plane(flex_points)
  n = matrix.col(p.normal)
  d = p.distance_to_origin
  assert approx_equal(abs(n), 1)
  dist0 = distance(n,d)
  for i in xrange(5000):
    d1 = d + random.uniform(-0.1, 0.1)
    n1 = matrix.rec(flex.random_double_r3_rotation_matrix(), (3,3))*n
    dist = distance(n1, d1)
    assert dist > dist0

def exercise_continued_fraction():
  continued_fraction = scitbx.math.continued_fraction
  frac = boost.rational.int
  cf = continued_fraction(1)
  assert cf.as_rational() == frac(1,1)
  cf.append(1)
  assert cf.as_rational() == frac(2,1)
  cf.append(1)
  assert cf.as_rational() == frac(3,2)
  cf.append(1)
  assert cf.as_rational() == frac(5,3)
  cf = continued_fraction.from_real(math.pi, eps=1e-2)
  assert cf.as_rational() == frac(22, 7)
  cf = continued_fraction.from_real(math.pi, eps=1e-4)
  assert cf.as_rational() == frac(333, 106)
  cf = continued_fraction.from_real(math.pi, eps=1e-5)
  assert cf.as_rational() == frac(355, 113)
  cf = continued_fraction.from_real(-math.pi, eps=1e-2)
  assert cf.as_rational() == frac(-22, 7)
  cf = continued_fraction.from_real(-math.pi, eps=1e-4)
  assert cf.as_rational() == frac(-333, 106)
  cf = continued_fraction.from_real(-math.pi, eps=1e-5)
  assert cf.as_rational() == frac(-355, 113)
  cf = continued_fraction.from_real(0.125)
  assert cf.as_rational() == frac(1,8)

def exercise_numeric_limits():
  l = scitbx.math.double_numeric_limits
  print "Floating point type 'double':"
  print "\tradix: ", l.radix
  print "\tmantissa digits (base 2):", l.digits
  print "\tmantissa digits (base 10):", l.digits10
  print "\tmin exponent (base 2):", l.min_exponent
  print "\tmin exponent (base 10):", l.min_exponent10
  print "\tmax exponent (base 2):", l.max_exponent
  print "\tmax exponent (base 10):", l.max_exponent10
  print "\tmin:", l.min
  print "\tmax:", l.max
  print "\tepsilon:", l.epsilon
  print "\tsafe min:", l.safe_min

def exercise_distributions():
  # normal distribution
  norm = distributions.normal_distribution()
  assert norm.mean() == 0
  assert norm.median() == 0
  assert norm.mode() == 0
  assert norm.standard_deviation() == 1
  assert norm.variance() == math.pow(norm.standard_deviation(), 2)
  assert norm.kurtosis() == 3
  assert norm.skewness() == 0
  assert approx_equal(norm.pdf(1.2), 0.19418605498321298)
  assert approx_equal(norm.cdf(norm.quantile(.9)), .9)
  assert approx_equal(norm.quantiles(5),
    (-1.2815515655446006, -0.52440051270804089, 0.0,
     0.52440051270804067, 1.2815515655446006))
  norm = distributions.normal_distribution(1,6)
  assert norm.mean() == 1
  assert norm.standard_deviation() == 6
  # student's t distribution
  try:
    stu = distributions.students_t_distribution(10)
  except RuntimeError, e:
    print "Skipping exercise students_t_distribution:", e
  else:
    assert stu.degrees_of_freedom() == 10
    assert stu.mean() == 0
    assert stu.median() == 0
    assert stu.mode() == 0
    assert approx_equal(
      math.pow(stu.standard_deviation(),2), 1.25)
    assert approx_equal(stu.variance(), 1.25)
    assert approx_equal(stu.kurtosis(), 4.0)
    assert stu.skewness() == 0
    assert approx_equal(norm.pdf(0.4), 0.066158757912835292)
    assert approx_equal(norm.cdf(norm.quantile(.8)), .8)
    assert approx_equal(stu.quantiles(6),
      (-1.4915762442496054, -0.69981206131243145, -0.21599563333226371,
       0.21599563333226388, 0.69981206131243145, 1.4915762442496057))

def exercise_approx_equal():
  from scitbx.math import double_numeric_limits as limits
  from scitbx.math import approx_equal_relatively

  # This would fail with a naive relative test for such tiny numbers
  assert approx_equal_relatively(-limits.min/2, limits.min/2,
                                 relative_error=1)
  # vanilla relative difference test
  assert approx_equal_relatively(0.9999, 1., 0.0001)
  assert approx_equal_relatively(0.9997 + 0.0004j, 1., 0.0005)

def exercise_weighted_covariance():
  from scitbx.math import weighted_covariance
  stats = weighted_covariance(x=flex.double((1, 2, 1e-4, 3, 4, 5)),
                              y=flex.double((2, 4, 1e-4, 6, 8, 10)),
                              weights=flex.double((2, 1, 3, 2, 1, 2)))
  eps = 1e-18
  # tests generated with Mathematica
  assert approx_equal(stats.mean_x, 2.18184545454545455, eps)
  assert approx_equal(stats.mean_y, 4.36366363636363636, eps)
  assert approx_equal(stats.variance_x, 3.42136859702479339, eps)
  assert approx_equal(stats.variance_y, 13.6857123986776860, eps)
  assert approx_equal(stats.covariance_xy, 6.84279669619834711, eps)
  assert approx_equal(stats.correlation, 0.99999999996534161, eps)
  stats.accumulate(x=-1, y=2, weight=2)
  stats.accumulate(x=2, y=-3, weight=1)
  assert approx_equal(stats.mean_x, 1.71430714285714286, eps)
  assert approx_equal(stats.mean_y, 3.50002142857142857, eps)
  assert approx_equal(stats.variance_x, 3.91829387923469388, eps)
  assert approx_equal(stats.variance_y, 14.6784214302551020, eps)
  assert approx_equal(stats.covariance_xy, 6.14274540984693878, eps)
  assert approx_equal(stats.correlation, 0.809980077408506317, eps)

def exercise_interpolation () :
  from scitbx.math import interpolate_catmull_rom_spline
  p0 = (4.6125, 53.1915, -1.0)
  p1 = (4.86, 54.206, 0.603)
  p2 = (6.640, 55.369, 0.651)
  p3 = (7.726, 56.192, -0.941)
  points = interpolate_catmull_rom_spline(p0, p1, p2, p3, 5)
  assert approx_equal(points[2][0], 5.9044, eps=0.0001)
  assert approx_equal(points[4][2], 0.651, eps=0.0001)
  p0 = (0, 1)
  p1 = (0.5, 0.707)
  p2 = (1, 0)
  p3 = (-0.5, -0.707)
  points = interpolate_catmull_rom_spline(p0, p1, p2, p3, 5)
  assert approx_equal(points[1][1], 0.454, eps=0.0001)

def run():
  exercise_weighted_covariance()
  exercise_distributions()
  exercise_approx_equal()
  exercise_median()
  exercise_numeric_limits()
  exercise_continued_fraction()
  exercise_least_squares_plane()
  exercise_div_mod()
  exercise_row_echelon_full_pivoting()
  exercise_solve_a_x_eq_b_min_norm_given_a_sym_b_col()
  exercise_eix()
  exercise_floating_point_epsilon()
  exercise_line_given_points()
  exercise_dihedral_angle()
  exercise_euler_angles()
  exercise_erf()
  exercise_gamma_incomplete()
  exercise_gamma_complete()
  exercise_exponential_integral_e1z()
  exercise_bessel()
  exercise_lambertw()
  exercise_golay()
  exercise_inertia_tensor()
  exercise_principal_axes_of_inertia()
  exercise_principal_axes_of_inertia_2d()
  explore_inertia_tensor_properties()
  exercise_phase_error()
  exercise_row_echelon()
  exercise_tensor_rank_2()
  exercise_icosahedron()
  exercise_basic_statistics()
  exercise_cheb_family()
  exercise_slatec_dlngam()
  exercise_slatec_dbinom()
  exercise_interpolation()
  forever = "--forever" in sys.argv[1:]
  exercise_unimodular_generator(
    forever=forever and "--unimodular" in sys.argv[1:])
  while 1:
    exercise_minimum_covering_sphere()
    if (not forever): break
  print "OK"

if (__name__ == "__main__"):
  run()
