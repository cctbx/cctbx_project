from scitbx.math import erf_verification, erf, erfc, erfcx
from scitbx.math import eigensystem
from scitbx.math import gaussian
from scitbx.array_family import flex
from scitbx.test_utils import approx_equal, eps_eq
import pickle
import StringIO
import random
import math

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
  g = gaussian.sum((1,-2,3,-4,5), (-.1,.2,-.3,.4,-.5))
  s = StringIO.StringIO()
  g.show(s)
  assert len(s.getvalue().split()) == 12

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
  gf = gaussian.fit(
    x, y, flex.double(),
    gaussian.sum((1,2), (4,5)))
  assert approx_equal(gf.table_sigmas(), [1,1,1])
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
  gf = gaussian.fit(
    x, reference_gaussian, flex.double(),
    gaussian.sum((1,2), (4,5)))
  assert approx_equal(gf.table_sigmas(), [1,1,1])
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
    x, reference_gaussian, flex.double(),
    gaussian.sum((1,2), (4,5)))
  sgf = gf.apply_shifts(flex.double((3,-3,4,6)), 0001)
  assert approx_equal(sgf.array_of_a(), (1+3,2+4))
  assert approx_equal(sgf.array_of_b(),
    ((math.sqrt(4)-3)**2,(math.sqrt(5)+6)**2))
  assert approx_equal(sgf.c(), 0)
  assert not sgf.use_c()
  sgf = gf.apply_shifts(flex.double((3,-3,4,6)), 00000)
  assert approx_equal(sgf.array_of_a(), (1+3,2+4))
  assert approx_equal(sgf.array_of_b(), (4-3,5+6))
  assert approx_equal(sgf.c(), 0)
  assert not sgf.use_c()
  differences = sgf.differences()
  for use_sigmas in [00000, 0001]:
    assert approx_equal(sgf.target_function(2, use_sigmas, differences),
      25.0320634)
    assert approx_equal(sgf.target_function(4, use_sigmas, differences),
      256.2682575)
    assert approx_equal(
      sgf.gradients_d_abc(2, use_sigmas, differences),
      [15.6539271, -4.1090114, 10.4562306, -1.6376781])
  gfc = gaussian.fit(
    x, reference_gaussian, flex.double(),
    gaussian.sum((1,2), (4,5), 6))
  sgfc = gfc.apply_shifts(flex.double((3,-3,4,6,-5)), 0001)
  assert approx_equal(sgfc.array_of_a(), (1+3,2+4))
  assert approx_equal(sgfc.array_of_b(),
    ((math.sqrt(4)-3)**2,(math.sqrt(5)+6)**2))
  assert approx_equal(sgfc.c(), 6-5)
  assert sgfc.use_c()
  sgfc = gfc.apply_shifts(flex.double((3,-3,4,6,-5)), 00000)
  assert approx_equal(sgfc.array_of_a(), (1+3,2+4))
  assert approx_equal(sgfc.array_of_b(), (4-3,5+6))
  assert approx_equal(sgfc.c(), 6-5)
  assert sgfc.use_c()
  differences = sgfc.differences()
  for use_sigmas in [00000, 0001]:
    assert approx_equal(sgfc.target_function(2, use_sigmas, differences),
      44.8181444)
    assert approx_equal(sgfc.target_function(4, use_sigmas, differences),
      757.3160329)
    assert approx_equal(
      sgfc.gradients_d_abc(2, use_sigmas, differences),
      [21.1132071, -6.0532695, 13.6638274, -2.2460994, 22.7860809])
  differences = c_fit.differences()
  gabc = c_fit.gradients_d_abc(2, 00000, differences)
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
  for include_constant_term in (00000, 0001):
    a = flex.double(g5c.array_of_a())
    b = flex.double(g5c.array_of_b())
    permutation = flex.sort_permutation(flex.abs(a), 1)[:4]
    gf = gaussian.fit(
      flex.double([0]),
      g5c,
      flex.double(),
      gaussian.sum(
        iter(a.select(permutation)),
        iter(b.select(permutation)), 0, include_constant_term))
    assert approx_equal(gf.differences(), [-5.01177418232])
    shifts = flex.double(8,-1)
    if (include_constant_term): shifts.append(-.2)
    sgf = gf.apply_shifts(shifts, 00000)
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
        flex.double(),
        sgf)
      differences = flex.double([0.5])
      assert approx_equal(
        gf.gradients_d_abc(2, 00000, differences),
        gaussian_fit_finite_diff_gradients(gf, gf.table_x()[0]),
        eps=1.e-3)
      for sigma in [0.04,0.02,0.01]:
        gf = gaussian.fit(
          flex.double([i / 20.]),
          g5c,
          flex.double([sigma]),
          sgf)
        for power in [2,4]:
          for use_sigmas in [00000, 0001]:
            differences = gf.differences()
            an=gf.gradients_d_abc(power, use_sigmas, differences)
            fi=gaussian_fit_finite_diff_target_gradients(gf, power, use_sigmas)
            assert eps_eq(an, fi, eps=1.e-3)

def run():
  exercise_erf()
  exercise_eigensystem()
  exercise_gaussian_term()
  exercise_gaussian_sum()
  exercise_gaussian_fit()
  print "OK"

if (__name__ == "__main__"):
  run()
