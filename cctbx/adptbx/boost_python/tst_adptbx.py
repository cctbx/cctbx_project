from cctbx import uctbx
from cctbx import adptbx
from cctbx import matrix
from scitbx.test_utils import approx_equal
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
  assert not adptbx.is_positive_definite(u)
  up = (0.534, 0.812, 0.613, 0.0166, 0.134, -0.0124)
  s = adptbx.eigensystem(up)
  s = adptbx.eigensystem(up, 1.e-7)
  assert approx_equal(s.values(), (0.813132, 0.432668, 0.713201))
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

def exercise_debye_waller():
  ucell = uctbx.unit_cell((5,7,9,80,100,130))
  assert abs(adptbx.debye_waller_factor_u_iso(ucell.d_star_sq((7,8,9))/4,0.025)
             - 0.006625427) < 1.e-6
  assert abs(adptbx.debye_waller_factor_u_iso(ucell, (7, 8, 9), 0.025)
             - 0.006625427) < 1.e-6
  u = (0.1, 0.2, 0.3, -0.04, 0.05, -0.06)
  assert abs(adptbx.debye_waller_factor_u_cif(ucell, (1, 2, 3), u)
             - 0.335822877927) < 1.e-6

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
    for i in xrange(3):
      check_eigenvector(u, es.values()[i], es.vectors(i))

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
  exercise_eigen()
  print "OK"

if (__name__ == "__main__"):
  run()
