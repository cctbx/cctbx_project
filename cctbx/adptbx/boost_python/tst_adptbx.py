from cctbx import uctbx
from cctbx import adptbx
from cctbx.array_family import flex
import scitbx.math.eigensystem
from scitbx import matrix
from libtbx.test_utils import approx_equal
import math
import random

def check_eigenvalue(adp, x):
  r = -adp[0] - adp[1] - adp[2]
  s =   adp[0] * adp[1] + adp[0] * adp[2] + adp[1] * adp[2] \
      - adp[3] * adp[3] - adp[4] * adp[4] - adp[5] * adp[5]
  t =   adp[0] * adp[5] * adp[5] - adp[0] * adp[1] * adp[2] \
      + adp[2] * adp[3] * adp[3] + adp[1] * adp[4] * adp[4] \
      - 2. * adp[3] * adp[4] * adp[5]
  f = x**3 + r * x**2 + s * x + t
  assert abs(f) < 1.e-4, f

def check_eigenvector(adp, x, v):
  v = matrix.col(v)
  assert approx_equal((matrix.sym(adp) * v * (1./x)).elems, v.elems, 1.e-4)

def exercise_interface():
  episq = 8*(math.pi**2)
  assert approx_equal(adptbx.u_as_b(2.3), 2.3*episq)
  assert approx_equal(adptbx.b_as_u(adptbx.u_as_b(2.3)), 2.3)
  u = (3,4,9, 2,1,7)
  assert approx_equal(adptbx.u_as_b(u), [x*episq for x in u])
  assert approx_equal(adptbx.b_as_u(adptbx.u_as_b(u)), u)
  uc = uctbx.unit_cell((5,4,7,80,110,100))
  for fw,bw in ((adptbx.u_cif_as_u_star, adptbx.u_star_as_u_cif),
                (adptbx.u_cart_as_u_star, adptbx.u_star_as_u_cart),
                (adptbx.u_cart_as_u_cif, adptbx.u_cif_as_u_cart),
                (adptbx.u_cart_as_beta, adptbx.beta_as_u_cart),
                (adptbx.u_cif_as_beta, adptbx.beta_as_u_cif)):
    assert approx_equal(bw(uc, fw(uc, u)), u)
  assert approx_equal(adptbx.beta_as_u_star(adptbx.u_star_as_beta(u)), u)
  assert approx_equal(adptbx.u_cart_as_u_iso(adptbx.u_iso_as_u_cart(2.3)), 2.3)
  for fw,bw in ((adptbx.u_iso_as_u_star, adptbx.u_star_as_u_iso),
                (adptbx.u_iso_as_u_cif, adptbx.u_cif_as_u_iso),
                (adptbx.u_iso_as_beta, adptbx.beta_as_u_iso)):
    assert approx_equal(bw(uc, fw(uc, 2.3)), 2.3)
  assert approx_equal(adptbx.debye_waller_factor_b_iso(0.25,2.3),
                      math.exp(-2.3*0.25))
  assert approx_equal(adptbx.debye_waller_factor_u_iso(0.25,2.3),
                      math.exp(-2.3*episq*0.25))
  assert approx_equal(adptbx.debye_waller_factor_b_iso(uc, (1,2,3), 2.3),
                      adptbx.debye_waller_factor_u_iso(uc, (1,2,3), 2.3/episq))
  u_star = adptbx.u_cart_as_u_star(uc, u)
  dw = adptbx.debye_waller_factor_u_star((1,2,3), u_star)
  assert approx_equal(dw, adptbx.debye_waller_factor_beta((1,2,3),
                            adptbx.u_star_as_beta(u_star)))
  assert approx_equal(dw, adptbx.debye_waller_factor_u_cif(uc, (1,2,3),
                            adptbx.u_star_as_u_cif(uc, u_star)))
  assert approx_equal(dw, adptbx.debye_waller_factor_u_cart(uc, (1,2,3),
                            adptbx.u_star_as_u_cart(uc, u_star)))
  for e in adptbx.eigenvalues(u):
    check_eigenvalue(u, e)
  assert not adptbx.is_positive_definite(adptbx.eigenvalues(u))
  assert not adptbx.is_positive_definite(adptbx.eigenvalues(u), 0)
  assert adptbx.is_positive_definite(adptbx.eigenvalues(u), 1.22)
  assert not adptbx.is_positive_definite(u)
  assert not adptbx.is_positive_definite(u, 0)
  assert adptbx.is_positive_definite(u, 1.22)
  up = (0.534, 0.812, 0.613, 0.0166, 0.134, -0.0124)
  s = adptbx.eigensystem(up)
  s = adptbx.eigensystem(up, 1.e-7)
  assert approx_equal(s.values(), (0.813132, 0.713201, 0.432668))
  for i in xrange(3):
    check_eigenvector(up, s.values()[i], s.vectors(i))
  c = (1,2,3, 3,-4,5, 4,5,6)
  v = (198,18,1020,116,447,269)
  assert approx_equal(adptbx.c_u_c_transpose(c, u), v)
  try: eigensystem(u)
  except: pass
  else: raise AssertionError, "Exception expected."
  s = adptbx.eigensystem(up)
  try: s.vectors(4)
  except: pass
  else: raise AssertionError, "Exception expected."
  uf = adptbx.eigenvalue_filtering(u)
  assert approx_equal(uf, (3.0810418, 4.7950710, 9.3400030,
                           1.7461615, 1.1659954, 6.4800706))
  assert approx_equal(scitbx.math.eigensystem.real_symmetric(u).values(),
                      (14.2792015, 2.9369144, -1.2161159))
  assert approx_equal(scitbx.math.eigensystem.real_symmetric(uf).values(),
                      (14.2792015, 2.9369144, 0))
  uf = adptbx.eigenvalue_filtering(up)
  assert approx_equal(uf, up)

def exercise_debye_waller():
  ucell = uctbx.unit_cell((5,7,9,80,100,130))
  assert abs(adptbx.debye_waller_factor_u_iso(ucell.d_star_sq((7,8,9))/4,0.025)
             - 0.006625427) < 1.e-6
  assert abs(adptbx.debye_waller_factor_u_iso(ucell, (7, 8, 9), 0.025)
             - 0.006625427) < 1.e-6
  u = (0.1, 0.2, 0.3, -0.04, 0.05, -0.06)
  assert abs(adptbx.debye_waller_factor_u_cif(ucell, (1, 2, 3), u)
             - 0.335822877927) < 1.e-6

def exercise_grad_u_transformations():
  uc = uctbx.unit_cell(
    (12.296512479074615, 15.985466222796999, 20.904071214426843,
     83.0, 109.0, 129.0))
  grad_u_star = (
    1681615.347859645, 1740095.5140965185, 2212142.6517873215,
    -2802250.0492254314, -1934508.748421974, 1503834.0901063806)
  grad_u_cart = adptbx.grad_u_star_as_u_cart(uc, grad_u_star)
  assert approx_equal(grad_u_cart,
    (11121.484290152286, 3713.1538936460488, 4293.4422553333607,
     -332.12536958297818, -340.3709419122527, 406.25522344235526))
  assert approx_equal(
    grad_u_star,
    adptbx.grad_u_cart_as_u_star(uc, grad_u_cart))
  grad_u_star_array = flex.sym_mat3_double([grad_u_star]*3)
  grad_u_cart_array = adptbx.grad_u_star_as_u_cart(uc, grad_u_star_array)
  for g in grad_u_cart_array:
    assert approx_equal(grad_u_cart, g)
  grad_u_star_array = adptbx.grad_u_cart_as_u_star(uc, grad_u_cart_array)
  for g in grad_u_star_array:
    assert approx_equal(grad_u_star, g)

def exercise_eigen_core(diag):
  u = adptbx.random_rotate_ellipsoid(diag + [0.,0.,0.])
  ev = list(adptbx.eigenvalues(u))
  diag.sort()
  ev.sort()
  for i in xrange(3):
    check_eigenvalue(u, ev[i])
  for i in xrange(3):
    assert abs(diag[i] - ev[i]) < 1.e-4
  if (adptbx.is_positive_definite(ev)):
    es = adptbx.eigensystem(u)
    ev = list(es.values())
    ev.sort()
    for i in xrange(3):
      check_eigenvalue(u, ev[i])
    for i in xrange(3):
      assert abs(diag[i] - ev[i]) < 1.e-4
    evec = []
    for i in xrange(3):
      check_eigenvector(u, es.values()[i], es.vectors(i))
      evec.extend(es.vectors(i))
    return # XXX following tests disabled for the moment
           # sometimes fail if eigenvalues are very similar but not identical
    sqrt_eval = matrix.diag(flex.sqrt(flex.double(es.values())))
    evec = matrix.sqr(evec).transpose()
    sqrt_u = evec * sqrt_eval * evec.transpose()
    u_full = matrix.sym(u).elems
    assert approx_equal(u_full, (sqrt_u.transpose()*sqrt_u).elems, eps=1.e-3)
    assert approx_equal(u_full, (sqrt_u*sqrt_u.transpose()).elems, eps=1.e-3)
    assert approx_equal(u_full, (sqrt_u*sqrt_u).elems, eps=1.e-3)
    sqrt_u_plus_shifts = matrix.sym(
      [x + 10*(random.random()-.5) for x in sqrt_u.as_sym_mat3()])
    sts = (sqrt_u_plus_shifts.transpose()*sqrt_u_plus_shifts).as_sym_mat3()
    ev = adptbx.eigenvalues(sts)
    assert min(ev) >= 0
    sts = (sqrt_u_plus_shifts*sqrt_u_plus_shifts.transpose()).as_sym_mat3()
    ev = adptbx.eigenvalues(sts)
    assert min(ev) >= 0
    sts = (sqrt_u_plus_shifts*sqrt_u_plus_shifts).as_sym_mat3()
    ev = adptbx.eigenvalues(sts)
    assert min(ev) >= 0

def exercise_eigen(n_trials=100):
  exercise_eigen_core([0,0,0])
  for n_equal in xrange(3):
    for i_trial in xrange(n_trials):
      if (n_equal == 0):
        diag = [random.random() for i in xrange(3)]
      if (n_equal == 1):
        i = random.randrange(3)
        diag[i] = random.random()
        diag[(i+1)%3] = random.random()
        diag[(i+2)%3] = diag[(i+1)%3]
      if (n_equal == 2):
        diag = [random.random()] * 3
      exercise_eigen_core(diag)

def run():
  exercise_interface()
  exercise_debye_waller()
  exercise_grad_u_transformations()
  exercise_eigen()
  print "OK"

if (__name__ == "__main__"):
  run()
