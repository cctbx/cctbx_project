import scitbx.math
from scitbx.math import euler_angles_as_matrix
from scitbx.math import erf_verification, erf, erfc, erfcx
from scitbx.math import bessel_i1_over_i0,bessel_i0,bessel_i1,bessel_ln_of_i0
from scitbx.math import bessel_inverse_i1_over_i0
from scitbx.math import matrix_inversion_in_place
from scitbx.math import eigensystem, time_eigensystem_real_symmetric
from scitbx.math import gaussian
from scitbx.math import golay_24_12_generator
from scitbx.math import principal_axes_of_inertia
from scitbx.math import sphere_3d, minimum_covering_sphere
from scitbx.math import signed_phase_error, phase_error, nearest_phase
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal, eps_eq
import pickle
import StringIO
import random
import math
import time
import sys

def exercise_floating_point_epsilon():
  float_eps = scitbx.math.floating_point_epsilon_float_get()
  double_eps = scitbx.math.floating_point_epsilon_double_get()
  assert 1.+float_eps != 1.
  assert 1.+double_eps != 1.
  assert float_eps >= double_eps
  assert 1.+double_eps/2. == 1.

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

def exercise_matrix_inversion_in_place():
  m = flex.double()
  m.resize(flex.grid(0,0))
  matrix_inversion_in_place(a=m)
  b = flex.double()
  b.resize(flex.grid(0,0))
  matrix_inversion_in_place(a=m, b=b)
  m = flex.double([2])
  m.resize(flex.grid(1,1))
  matrix_inversion_in_place(m)
  assert approx_equal(m, [1/2.])
  m = flex.double([2,0,0,-3])
  m.resize(flex.grid(2,2))
  matrix_inversion_in_place(m)
  assert approx_equal(m, [1/2.,0,0,-1/3.])
  m = flex.double([1,2,-3,4])
  m.resize(flex.grid(2,2))
  matrix_inversion_in_place(m)
  assert approx_equal(m, [2/5.,-1/5.,3/10.,1/10.])
  m = flex.double([2,0,0,0,-3,0,0,0,4])
  m.resize(flex.grid(3,3))
  from scitbx import matrix
  matrix_inversion_in_place(m)
  assert approx_equal(m, [1/2.,0,0,0,-1/3.,0,0,0,1/4.])
  m = flex.double([1,2,-3,-2,4,-1,8,0,4])
  m.resize(flex.grid(3,3))
  matrix_inversion_in_place(m)
  assert approx_equal(m, [1/7.,-1/14.,5/56.,0,1/4.,1/16.,-2/7.,1/7.,1/14.])
  for n in xrange(1,12):
    u = flex.double(n*n, 0)
    for i in xrange(0,n*n,n+1): u[i] = 1
    for diag in [1,2]:
      m = flex.double(n*n, 0)
      for i in xrange(0,n*n,n+1): m[i] = diag
      m.resize(flex.grid(n,n))
      m_orig = matrix.rec(m, (n,n))
      matrix_inversion_in_place(m)
      m_inv = matrix.rec(m, (n,n))
      assert approx_equal(m_orig*m_inv, u)
      assert approx_equal(m_inv*m_orig, u)
      for n_b in xrange(0,4):
        m = flex.double(m_orig)
        m.resize(flex.grid(n,n))
        b = flex.double(xrange(1,n*n_b+1))
        b.resize(flex.grid(n_b,n))
        b_orig = matrix.rec(b, (n,n_b))
        matrix_inversion_in_place(m, b)
        m_inv = matrix.rec(m, (n,n))
        x = matrix.rec(b, (n_b,n))
        assert approx_equal(m_orig*m_inv, u)
        assert approx_equal(m_inv*m_orig, u)
        for i_b in xrange(n_b):
          b_i = matrix.col(b_orig.elems[i_b*n:(i_b+1)*n])
          x_i = matrix.col(x.elems[i_b*n:(i_b+1)*n])
          assert approx_equal(m_orig*x_i, b_i)
  for n in xrange(1,12):
    u = flex.double(n*n, 0)
    for i in xrange(0,n*n,n+1): u[i] = 1
    for i_trial in xrange(3):
      m = 2*flex.random_double(n*n)-1
      m.resize(flex.grid(n,n))
      m_orig = matrix.rec(m, (n,n))
      try:
        matrix_inversion_in_place(m)
      except RuntimeError, e:
        assert str(e) == "scitbx Error: matrix is singular."
      else:
        m_inv = matrix.rec(m, (n,n))
        assert approx_equal(m_orig*m_inv, u)
        assert approx_equal(m_inv*m_orig, u)
        for n_b in xrange(0,4):
          m = flex.double(m_orig)
          m.resize(flex.grid(n,n))
          b = flex.random_double(n*n_b)
          b.resize(flex.grid(n_b,n))
          b_orig = matrix.rec(b, (n,n_b))
          matrix_inversion_in_place(m, b)
          m_inv = matrix.rec(m, (n,n))
          x = matrix.rec(b, (n_b,n))
          assert approx_equal(m_orig*m_inv, u)
          assert approx_equal(m_inv*m_orig, u)
          for i_b in xrange(n_b):
            b_i = matrix.col(b_orig.elems[i_b*n:(i_b+1)*n])
            x_i = matrix.col(x.elems[i_b*n:(i_b+1)*n])
            assert approx_equal(m_orig*x_i, b_i)

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

def exercise_eigensystem():
  for n in xrange(1,10):
    m = flex.double(flex.grid(n,n))
    s = eigensystem.real_symmetric(m)
    assert approx_equal(tuple(s.values()), [0]*n)
    v = s.vectors()
    for i in xrange(n):
      for j in xrange(n):
        x = 0
        if (i == j): x = 1
        assert approx_equal(v[(i,j)], x)
    v = []
    for i in xrange(n):
      j = (i*13+17) % n
      v.append(j)
      m[i*(n+1)] = j
    s = eigensystem.real_symmetric(m)
    if (n == 3):
      ss = eigensystem.real_symmetric((m[0],m[4],m[8],m[1],m[2],m[5]))
      assert approx_equal(s.values(), ss.values())
      assert approx_equal(s.vectors(), ss.vectors())
    v.sort()
    v.reverse()
    assert approx_equal(s.values(), v)
    if (n > 1):
      assert approx_equal(flex.min(s.vectors()), 0)
    assert approx_equal(flex.max(s.vectors()), 1)
    assert approx_equal(flex.sum(s.vectors()), n)
    for t in xrange(10):
      for i in xrange(n):
        for j in xrange(i,n):
          m[i*n+j] = random.random() - 0.5
          if (i != j):
            m[j*n+i] = m[i*n+j]
      s = eigensystem.real_symmetric(m)
      if (n == 3):
        ss = eigensystem.real_symmetric((m[0],m[4],m[8],m[1],m[2],m[5]))
        assert approx_equal(s.values(), ss.values())
        assert approx_equal(s.vectors(), ss.vectors())
      v = list(s.values())
      v.sort()
      v.reverse()
      assert list(s.values()) == v
      for i in xrange(n):
        l = s.values()[i]
        x = s.vectors()[i*n:i*n+n]
        mx = matrix_mul(m, n, n, x, n, 1)
        lx = [e*l for e in x]
        assert approx_equal(mx, lx)
  m = (1.4573362052597449, 1.7361052947659894, 2.8065584999742659,
       -0.5387293498219814, -0.018204949672480729, 0.44956507395617257)
  n_repetitions = 100000
  t0 = time.time()
  v = time_eigensystem_real_symmetric(m, n_repetitions)
  assert v == (0,0,0)
  print "time_eigensystem_real_symmetric: %.3f micro seconds" % (
    (time.time() - t0)/n_repetitions*1.e6)

def gaussian_finite_gradient_dx_at_x(gaussian, x, eps=1.e-5):
  if (x == 0): return 0
  assert x >= eps
  tm = gaussian.at_x(x-eps)
  tp = gaussian.at_x(x+eps)
  return (tp-tm)/(2*eps)

def exercise_gaussian_gradient_dx(gaussian, x_max=1., n_points=50):
  for i in xrange(n_points+1):
    x = x_max * i / n_points
    grad_finite = gaussian_finite_gradient_dx_at_x(gaussian, x)
    grad_analytical = gaussian.gradient_dx_at_x(x)
    assert eps_eq(grad_finite, grad_analytical)

def exercise_gaussian_integral_dx(gaussian, x_max=1., n_points=1000):
  numerical_integral = 0
  x_step = x_max / n_points
  for i in xrange(n_points+1):
    x = x_max * i / n_points
    new_value = gaussian.at_x(x)
    if (i):
      numerical_integral += (prev_value + new_value) * .5
    prev_value = new_value
    analytical_integral = gaussian.integral_dx_at_x(x, 1.e-3)
    assert eps_eq(analytical_integral, gaussian.integral_dx_at_x(x))
    assert eps_eq(numerical_integral*x_step, analytical_integral, eps=1.e-5)

def gaussian_term_finite_gradient_d_ab_at_x(term, x, eps=1.e-5):
  tm = gaussian.term(term.a-eps,term.b).at_x(x)
  tp = gaussian.term(term.a+eps,term.b).at_x(x)
  gr_a = (tp-tm)/(2*eps)
  tm = gaussian.term(term.a,term.b-eps).at_x(x)
  tp = gaussian.term(term.a,term.b+eps).at_x(x)
  gr_b = (tp-tm)/(2*eps)
  return gaussian.term(gr_a, gr_b)

def exercise_gaussian_term_gradients_d_ab(term, x_max=1., n_points=50):
  for i in xrange(n_points+1):
    x = x_max * i / n_points
    grad_finite = gaussian_term_finite_gradient_d_ab_at_x(term, x)
    grad_analytical = term.gradients_d_ab_at_x_sq(x*x)
    assert eps_eq(grad_finite.a, grad_analytical.a)
    assert eps_eq(grad_finite.b, grad_analytical.b)

def exercise_gaussian_term():
  t = gaussian.term(2,3)
  assert approx_equal(t.a, 2)
  assert approx_equal(t.b, 3)
  assert approx_equal(t.at_x_sq(4), 2*math.exp(-3*4))
  assert approx_equal(t.at_x(2), 2*math.exp(-3*4))
  eps = 1.e-5
  for ix in (xrange(10)):
    x = ix/10.
    assert eps_eq((t.at_x(x+eps)-t.at_x(x-eps))/(2*eps), t.gradient_dx_at_x(x))
  for f in [1,-1]:
    for t in [gaussian.term(f*2,3),
              gaussian.term(f*3,0),
              gaussian.term(f*4,1.e-4),
              gaussian.term(f*5,-1)]:
      exercise_gaussian_gradient_dx(t)
      exercise_gaussian_integral_dx(t)
      exercise_gaussian_term_gradients_d_ab(t)

def exercise_gaussian_sum():
  g = gaussian.sum(0)
  assert g.n_terms() == 0
  assert g.array_of_a() == ()
  assert g.array_of_b() == ()
  assert approx_equal(g.c(), 0)
  assert g.use_c()
  assert g.n_parameters() == 1
  g = gaussian.sum(0, True)
  assert g.use_c()
  g = gaussian.sum(0, False)
  assert not g.use_c()
  g = gaussian.sum(1)
  assert g.n_terms() == 0
  assert g.array_of_a() == ()
  assert g.array_of_b() == ()
  assert approx_equal(g.c(), 1)
  assert g.use_c()
  assert g.n_parameters() == 1
  g = gaussian.sum((), ())
  assert g.n_terms() == 0
  assert g.array_of_a() == ()
  assert g.array_of_b() == ()
  assert g.c() == 0
  assert not g.use_c()
  assert g.n_parameters() == 0
  g = gaussian.sum((), (), -2)
  assert g.n_terms() == 0
  assert g.array_of_a() == ()
  assert g.array_of_b() == ()
  assert approx_equal(g.c(), -2)
  g = gaussian.sum(flex.double((1,2,3,4)))
  assert approx_equal(g.array_of_a(), (1,3))
  assert approx_equal(g.array_of_b(), (2,4))
  assert approx_equal(g.c(), 0)
  assert not g.use_c()
  g = gaussian.sum(flex.double((1,2,3,4)), 0, True)
  assert approx_equal(g.c(), 0)
  assert g.use_c()
  g = gaussian.sum(flex.double((1,2,3,4)), 5)
  assert approx_equal(g.c(), 5)
  assert g.use_c()
  g = gaussian.sum((1,-2,3,-4,5), (-.1,.2,-.3,.4,-.5), 6)
  assert g.n_terms() == 5
  assert approx_equal(g.array_of_a(),(1,-2,3,-4,5))
  assert approx_equal(g.array_of_b(),(-.1,.2,-.3,.4,-.5))
  assert approx_equal(g.c(), 6)
  assert approx_equal(g.at_x_sq(3/4.), 13.4251206)
  assert approx_equal(g.at_x_sq(flex.double([2/4.,3/4.])),
    [11.8723031, 13.4251206])
  assert approx_equal(g.at_x(math.sqrt(3/4.)), 13.4251206)
  assert approx_equal(g.at_x(flex.sqrt(flex.double([2/4.,3/4.]))),
    [11.8723031, 13.4251206])
  s = pickle.dumps(g)
  l = pickle.loads(s)
  assert l.n_terms() == g.n_terms()
  assert approx_equal(l.array_of_a(), g.array_of_a())
  assert approx_equal(l.array_of_b(), g.array_of_b())
  assert approx_equal(l.c(), g.c())
  assert l.use_c()
  s = pickle.dumps(gaussian.sum((),()))
  l = pickle.loads(s)
  assert not l.use_c()
  exercise_gaussian_gradient_dx(gaussian.sum(
    [5.5480], [10.4241], 0))
  exercise_gaussian_gradient_dx(gaussian.sum(
   [2.657506,1.078079,1.490909,-4.241070,0.713791],
   [14.780758,0.776775,42.086842,-0.000294,0.239535],
   4.297983))
  exercise_gaussian_integral_dx(gaussian.sum(
    [5.5480], [10.4241]))
  exercise_gaussian_integral_dx(gaussian.sum(
    [5.5480], [10.4241], 3))
  exercise_gaussian_integral_dx(gaussian.sum(
    [5.5480], [0], 0))
  exercise_gaussian_integral_dx(gaussian.sum(
    [5.5480], [-0.01]))
  exercise_gaussian_integral_dx(gaussian.sum(
    [2.657506,1.078079,1.490909,-4.241070,0.713791],
    [14.780758,0.776775,42.086842,-0.000294,0.239535],
    4.297983))
  g = gaussian.sum((1,-2,3,-4,5), (-.1,.2,-.3,.4,-.5), 6)
  s = StringIO.StringIO()
  g.show(s)
  assert len(s.getvalue().split()) == 14
  g = gaussian.sum((3,-2,1,-4,5), (-.3,.2,-.1,.4,-.5))
  s = StringIO.StringIO()
  g.show(s)
  assert len(s.getvalue().split()) == 12
  assert isinstance(g.sort(), gaussian.sum)
  assert approx_equal(g.sort().array_of_a(), (5,-4,3,-2,1))
  assert approx_equal(g.sort().array_of_b(), (-.5,.4,-.3,.2,-.1))
  assert not g.sort().use_c()
  g = gaussian.sum((1,2),(3,4),5)
  assert approx_equal(g.sort().array_of_a(), (2,1))
  assert approx_equal(g.sort().array_of_b(), (4,3))
  assert approx_equal(g.sort().c(), 5)
  assert g.sort().use_c()

def gaussian_fit_finite_diff_gradients(gfit, x, eps=1.e-2):
  gr = flex.double()
  c = gfit.c()
  use_c = gfit.use_c()
  for i in xrange(gfit.n_terms()):
    t = []
    for seps in (eps, -eps):
      a = list(gfit.array_of_a())
      a[i] += seps
      t.append(
        gaussian.sum(a, gfit.array_of_b(), c, use_c).at_x(x))
    gr.append((t[0]-t[1])/(2*eps))
    t = []
    for seps in (eps, -eps):
      b = list(gfit.array_of_b())
      b[i] += seps
      t.append(
        gaussian.sum(gfit.array_of_a(), b, c, use_c).at_x(x))
    gr.append((t[0]-t[1])/(2*eps))
  if (use_c):
    t = []
    for seps in (eps, -eps):
      t.append(
        gaussian.sum(
          gfit.array_of_a(), gfit.array_of_b(), c+seps, use_c).at_x(x))
    gr.append((t[0]-t[1])/(2*eps))
  return gr

def gaussian_fit_finite_diff_target_gradients(gfit, power, use_sigmas,
                                              eps=1.e-2):
  assert gfit.table_x().size() == 1
  weight = 1/gfit.table_sigmas()[0]**2
  gr = flex.double()
  c = gfit.c()
  use_c = gfit.use_c()
  for i in xrange(gfit.n_terms()):
    t = []
    for seps in (eps, -eps):
      a = list(gfit.array_of_a())
      a[i] += seps
      gf = gaussian.fit(
        gfit.table_x(),
        gfit.table_y(),
        gfit.table_sigmas(),
        gaussian.sum(a, gfit.array_of_b(), c, use_c))
      t.append(gf.target_function(power, use_sigmas, gf.differences()))
    gr.append((t[0]-t[1])/(2*eps))
    t = []
    for seps in (eps, -eps):
      b = list(gfit.array_of_b())
      b[i] += seps
      gf = gaussian.fit(
        gfit.table_x(),
        gfit.table_y(),
        gfit.table_sigmas(),
        gaussian.sum(gfit.array_of_a(), b, c, use_c))
      t.append(gf.target_function(power, use_sigmas, gf.differences()))
    gr.append((t[0]-t[1])/(2*eps))
  if (use_c):
    t = []
    for seps in (eps, -eps):
      gf = gaussian.fit(
        gfit.table_x(),
        gfit.table_y(),
        gfit.table_sigmas(),
        gaussian.sum(gfit.array_of_a(), gfit.array_of_b(), c+seps, use_c))
      t.append(gf.target_function(power, use_sigmas, gf.differences()))
    gr.append((t[0]-t[1])/(2*eps))
  return gr

def exercise_gaussian_fit():
  x = flex.double((0.1, 0.2, 0.5))
  y = flex.double((3,2,1))
  sigmas = flex.double((0.04,0.02,0.01))
  gf = gaussian.fit(
    x, y, sigmas,
    gaussian.sum((1,2), (4,5)))
  assert approx_equal(gf.array_of_a(), (1,2))
  assert approx_equal(gf.array_of_b(), (4,5))
  assert approx_equal(gf.c(), 0)
  assert not gf.use_c()
  assert approx_equal(gf.table_x(), x)
  assert approx_equal(gf.table_y(), y)
  assert approx_equal(gf.table_sigmas(), sigmas)
  assert approx_equal(gf.fitted_values(),
    [2.8632482881537511, 2.4896052951221748, 0.94088903489182252])
  reference_gaussian = gaussian.sum((1,2,3), (4,5,6))
  gf = gaussian.fit(
    x, reference_gaussian, sigmas,
    gaussian.sum((1,2), (4,5)))
  assert approx_equal(gf.array_of_a(), (1,2))
  assert approx_equal(gf.array_of_b(), (4,5))
  assert approx_equal(gf.c(), 0)
  assert approx_equal(gf.table_x(), x)
  assert approx_equal(gf.table_y(), reference_gaussian.at_x(x))
  assert approx_equal(gf.table_sigmas(), sigmas)
  assert isinstance(gf.sort(), gaussian.fit)
  assert gf.sort().table_x() == gf.table_x()
  assert gf.sort().table_y() == gf.table_y()
  assert gf.sort().table_sigmas() == gf.table_sigmas()
  assert approx_equal(gf.differences(), gf.at_x(x)-reference_gaussian.at_x(x))
  c_fit = gaussian.fit(
    flex.double([0.0, 0.066666666666666666, 0.13333333333333333,
                 0.2, 0.26666666666666666]),
    gaussian.sum(
      (2.657506, 1.078079, 1.490909, -4.2410698, 0.71379101),
      (14.780758, 0.776775, 42.086842, -0.000294, 0.239535),
      4.2979832),
    flex.double(5, 0.0005),
    gaussian.sum(
      (1.1423916, 4.1728425, 0.61716694),
      (0.50733125, 14.002512, 41.978928)))
  differences = flex.double([-0.064797341823577881, 0.003608505180995536,
    0.098159179757290715, 0.060724224581695019, -0.10766283796372011])
  assert approx_equal(c_fit.differences(), differences)
  assert approx_equal(c_fit.significant_relative_errors(),
    [0.0107212, 0.0005581, 0.0213236, 0.0169304, 0.0385142])
  gf = gaussian.fit(
    x, reference_gaussian, flex.double(x.size(), 1),
    gaussian.sum((1,2), (4,5)))
  assert list(gf.bound_flags(False, False)) == [False,False,False,False]
  assert list(gf.bound_flags(True, False)) == [True,False,True,False]
  assert list(gf.bound_flags(False, True)) == [False,True,False,True]
  sgf = gf.apply_shifts(flex.double((3,-3,4,6)), True)
  assert approx_equal(sgf.array_of_a(), (1+3,2+4))
  assert approx_equal(sgf.array_of_b(),
    ((math.sqrt(4)-3)**2,(math.sqrt(5)+6)**2))
  assert approx_equal(sgf.c(), 0)
  assert not sgf.use_c()
  sgf = gf.apply_shifts(flex.double((3,-3,4,6)), False)
  assert approx_equal(sgf.array_of_a(), (1+3,2+4))
  assert approx_equal(sgf.array_of_b(), (4-3,5+6))
  assert approx_equal(sgf.c(), 0)
  assert not sgf.use_c()
  differences = sgf.differences()
  for use_sigmas in [False, True]:
    assert approx_equal(sgf.target_function(2, use_sigmas, differences),
      25.0320634)
    assert approx_equal(sgf.target_function(4, use_sigmas, differences),
      256.2682575)
    assert approx_equal(
      sgf.gradients_d_abc(2, use_sigmas, differences),
      [15.6539271, -4.1090114, 10.4562306, -1.6376781])
  gfc = gaussian.fit(
    x, reference_gaussian, flex.double(x.size(), 1),
    gaussian.sum((1,2), (4,5), 6))
  assert list(gfc.bound_flags(False, False)) == [False,False,False,False,False]
  assert list(gfc.bound_flags(True, False)) == [True,False,True,False,True]
  assert list(gfc.bound_flags(False, True)) == [False,True,False,True,False]
  sgfc = gfc.apply_shifts(flex.double((3,-3,4,6,-5)), True)
  assert approx_equal(sgfc.array_of_a(), (1+3,2+4))
  assert approx_equal(sgfc.array_of_b(),
    ((math.sqrt(4)-3)**2,(math.sqrt(5)+6)**2))
  assert approx_equal(sgfc.c(), 6-5)
  assert sgfc.use_c()
  sgfc = gfc.apply_shifts(flex.double((3,-3,4,6,-5)), False)
  assert approx_equal(sgfc.array_of_a(), (1+3,2+4))
  assert approx_equal(sgfc.array_of_b(), (4-3,5+6))
  assert approx_equal(sgfc.c(), 6-5)
  assert sgfc.use_c()
  differences = sgfc.differences()
  for use_sigmas in [False, True]:
    assert approx_equal(sgfc.target_function(2, use_sigmas, differences),
      44.8181444)
    assert approx_equal(sgfc.target_function(4, use_sigmas, differences),
      757.3160329)
    assert approx_equal(
      sgfc.gradients_d_abc(2, use_sigmas, differences),
      [21.1132071, -6.0532695, 13.6638274, -2.2460994, 22.7860809])
  differences = c_fit.differences()
  gabc = c_fit.gradients_d_abc(2, False, differences)
  assert approx_equal(
    gabc,
    [-0.016525391425206391, 0.0074465239375589107, 0.020055876723667564,
     0.00054794635257838251, -0.018754011379726425, -0.0011194004809549143])
  assert approx_equal(
    c_fit.gradients_d_shifts(flex.double((0.1,0.4,0.2,0.5,0.3,0.6)), gabc),
    [-0.0165254, 0.01656512, 0.0200559, 0.0046488, -0.0187540, -0.0158487])
  g5c = gaussian.sum(
    (2.657505989074707, 1.0780789852142334, 1.4909089803695679,
     -4.2410697937011719, 0.71379101276397705),
    (14.780757904052734, 0.77677500247955322, 42.086841583251953,
     -0.00029399999766610563, 0.23953500390052795),
    4.2979831695556641)
  for include_constant_term in (False, True):
    a = flex.double(g5c.array_of_a())
    b = flex.double(g5c.array_of_b())
    permutation = flex.sort_permutation(data=flex.abs(a), reverse=True)[:4]
    gf = gaussian.fit(
      flex.double([0]),
      g5c,
      flex.double(1, 1),
      gaussian.sum(
        iter(a.select(permutation)),
        iter(b.select(permutation)), 0, include_constant_term))
    assert approx_equal(gf.differences(), [-5.01177418232])
    shifts = flex.double(8,-1)
    if (include_constant_term): shifts.append(-.2)
    sgf = gf.apply_shifts(shifts, False)
    assert approx_equal(sgf.array_of_a(),
                        [-5.2410698, 1.657506, 0.49090898, 0.078078985])
    assert approx_equal(sgf.array_of_b(),
                        [-1.0002940, 13.780758, 41.086842, -0.223225])
    if (include_constant_term):
      assert approx_equal(sgf.c(), -.2)
    expected_gradients = [1,0,1,0,1,0,1,0]
    if (include_constant_term): expected_gradients.append(1)
    assert approx_equal(
      gaussian_fit_finite_diff_gradients(sgf, 0),
      expected_gradients,
      eps=1.e-4)
    for i in xrange(10):
      gf = gaussian.fit(
        flex.double([i / 10.]),
        g5c,
        flex.double(1, 1),
        sgf)
      differences = flex.double([0.5])
      assert approx_equal(
        gf.gradients_d_abc(2, False, differences),
        gaussian_fit_finite_diff_gradients(gf, gf.table_x()[0]),
        eps=1.e-3)
      for sigma in [0.04,0.02,0.01]:
        gf = gaussian.fit(
          flex.double([i / 20.]),
          g5c,
          flex.double([sigma]),
          sgf)
        for power in [2,4]:
          for use_sigmas in [False, True]:
            differences = gf.differences()
            an=gf.gradients_d_abc(power, use_sigmas, differences)
            fi=gaussian_fit_finite_diff_target_gradients(gf, power, use_sigmas)
            assert eps_eq(an, fi, eps=1.e-3)

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
    raise RuntimeError("Exception expected.")
  weights = [0]*25
  for code in golay_24_12_generator():
    weights[list(code).count(1)] += 1
  assert weights == [1,0,0,0,0,0,0,0,759,0,0,0,2576,0,0,0,759,0,0,0,0,0,0,0,1]

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
  for i_trial in xrange(10):
    if (i_trial == 0):
      center_of_mass = [0,0,0]
    else:
      center_of_mass = [rnd(),rnd(),rnd()]
    points = flex.vec3_double()
    for point in flex.nested_loop([-1,-1,-1], [2,2,2]):
      points.append((matrix.col(point) + matrix.col(center_of_mass)).elems)
    pai = principal_axes_of_inertia(points=points)
    assert approx_equal(pai.center_of_mass(), center_of_mass)
    assert approx_equal(pai.inertia_tensor(), [36,36,36,0,0,0])
    es = pai.eigensystem()
    assert approx_equal(es.values(), [36,36,36])
    if (i_trial == 0):
      assert approx_equal(es.vectors(), [1,0,0,0,1,0,0,0,1])
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

def exercise_minimum_covering_sphere(epsilon=1.e-3):
  s3 = sphere_3d(center=[1,2,3], radius=4)
  assert approx_equal(s3.center(), [1,2,3])
  assert approx_equal(s3.radius(), 4)
  s3 = s3.expand(additional_radius=2)
  assert approx_equal(s3.center(), [1,2,3])
  assert approx_equal(s3.radius(), 6)
  assert s3.is_inside(point=[1,2,3])
  assert s3.is_inside(point=[1,2,3+6-1.e-6])
  assert not s3.is_inside([1,2,3+6+1.e-6])
  assert approx_equal(s3.box_min(), [1-6,2-6,3-6])
  assert approx_equal(s3.box_max(), [1+6,2+6,3+6])
  points = flex.vec3_double([(0,0,0),(1,0,0),(0,1,0),(1,1,1)])
  mcs = scitbx.math.minimum_covering_sphere_3d(points=points)
  assert mcs.n_iterations() > 0
  assert approx_equal(mcs.center(), (0.5,0.5,0.5), eps=1.e-3)
  assert approx_equal(mcs.radius(), math.sqrt(3)/2, eps=1.e-5)
  assert mcs.is_inside(mcs.center()) # base class method
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

def run():
  exercise_floating_point_epsilon()
  exercise_euler_angles()
  exercise_erf()
  exercise_bessel()
  exercise_matrix_inversion_in_place()
  exercise_eigensystem()
  exercise_gaussian_term()
  exercise_gaussian_sum()
  exercise_gaussian_fit()
  exercise_golay()
  exercise_principal_axes_of_inertia()
  exercise_phase_error()
  forever = "--Forever" in sys.argv[1:]
  while 1:
    exercise_minimum_covering_sphere()
    if (not forever): break
  print "OK"

if (__name__ == "__main__"):
  run()
