from __future__ import absolute_import, division, print_function
from cctbx import uctbx
from cctbx import adptbx
from cctbx.array_family import flex
import scitbx.linalg.eigensystem
from scitbx import matrix
from libtbx.test_utils import Exception_expected, approx_equal
import math
import random
from six.moves import range
from six.moves import zip

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
  assert approx_equal(
    (matrix.sym(sym_mat3=adp) * v * (1./x)).elems, v.elems, 1.e-4)

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
  fc = adptbx.factor_u_cart_u_iso(u_cart=u)
  assert approx_equal(fc.u_iso, adptbx.u_cart_as_u_iso(u))
  assert approx_equal(
    fc.u_cart_minus_u_iso,
    [uii-fc.u_iso for uii in u[:3]]+list(u[3:]))
  f = adptbx.factor_u_star_u_iso(
    unit_cell=uc, u_star=adptbx.u_cart_as_u_star(uc, u))
  assert approx_equal(f.u_iso, fc.u_iso)
  assert approx_equal(
    f.u_star_minus_u_iso,
    adptbx.u_cart_as_u_star(uc, fc.u_cart_minus_u_iso))
  f = adptbx.factor_u_cif_u_iso(
    unit_cell=uc, u_cif=adptbx.u_cart_as_u_cif(uc, u))
  assert approx_equal(f.u_iso, fc.u_iso)
  assert approx_equal(
    f.u_cif_minus_u_iso,
    adptbx.u_cart_as_u_cif(uc, fc.u_cart_minus_u_iso))
  f = adptbx.factor_beta_u_iso(
    unit_cell=uc, beta=adptbx.u_cart_as_beta(uc, u))
  assert approx_equal(f.u_iso, fc.u_iso)
  assert approx_equal(
    f.beta_minus_u_iso,
    adptbx.u_cart_as_beta(uc, fc.u_cart_minus_u_iso))
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
  assert approx_equal(s.values(), (0.813132, 0.713201, 0.432668))
  for i in range(3):
    check_eigenvector(up, s.values()[i], s.vectors(i))
  c = (1,2,3, 3,-4,5, 4,5,6)
  v = (198,18,1020,116,447,269)
  assert approx_equal(adptbx.c_u_c_transpose(c, u), v)
  assert approx_equal(adptbx.eigensystem(u).values(),
    (14.279201519086316, 2.9369143826320214, -1.2161159017183376))
  s = adptbx.eigensystem(up)
  try: s.vectors(4)
  except RuntimeError as e: assert str(e).endswith("Index out of range.")
  else: raise Exception_expected
  uf = adptbx.eigenvalue_filtering(u_cart=u, u_min=0)
  assert approx_equal(uf, (3.0810418, 4.7950710, 9.3400030,
                           1.7461615, 1.1659954, 6.4800706))
  uf = adptbx.eigenvalue_filtering(u_cart=u, u_min=0, u_max=3)
  assert approx_equal(uf, (2.7430890, 1.0378360, 2.1559895,
                           0.6193215, -0.3921632, 1.2846854))
  uf = adptbx.eigenvalue_filtering(u_cart=u, u_min=0, u_max=3)
  assert approx_equal(scitbx.linalg.eigensystem.real_symmetric(u).values(),
                      (14.2792015, 2.9369144, -1.2161159))
  assert approx_equal(scitbx.linalg.eigensystem.real_symmetric(uf).values(),
                      (3, 2.9369144, 0))
  uf = adptbx.eigenvalue_filtering(up)
  assert approx_equal(uf, up)

def u_star_minus_u_iso_airlie(unit_cell, u_star):
  u_iso = adptbx.u_star_as_u_iso(unit_cell,u_star)
  u_star_as_beta = adptbx.u_star_as_beta(u_star)
  u_iso_as_beta = adptbx.u_iso_as_beta(unit_cell,u_iso)
  beta_minus_u_iso = [a-i for a,i in zip(u_star_as_beta,u_iso_as_beta)]
  return adptbx.beta_as_u_star(beta_minus_u_iso)

def u_star_minus_u_iso_ralf(unit_cell, u_star):
  u_cart = adptbx.u_star_as_u_cart(unit_cell, u_star)
  u_iso = adptbx.u_cart_as_u_iso(u_cart)
  u_cart_minus_u_iso = [a-u_iso for a in u_cart[:3]] + list(u_cart[3:])
  return adptbx.u_cart_as_u_star(unit_cell, u_cart_minus_u_iso)

def exercise_factor_u_star_u_iso():
  for i_trial in range(100):
    a = flex.random_double(size=9, factor=3)
    a.resize(flex.grid(3,3))
    u = a.matrix_transpose().matrix_multiply(a) # always positive-definite
    u_cart = [u[0],u[4],u[8],u[1],u[2],u[5]]
    unit_cell = uctbx.unit_cell((3,5,7,80,100,110))
    u_star = adptbx.u_cart_as_u_star(unit_cell, u_cart)
    airlie = u_star_minus_u_iso_airlie(unit_cell, u_star)
    ralf = u_star_minus_u_iso_ralf(unit_cell, u_star)
    assert approx_equal(ralf, airlie, 1.e-10)
    f = adptbx.factor_u_star_u_iso(unit_cell=unit_cell, u_star=u_star)
    assert approx_equal(f.u_iso, adptbx.u_cart_as_u_iso(u_cart))
    assert approx_equal(f.u_star_minus_u_iso, airlie, 1.e-10)

def exercise_debye_waller():
  ucell = uctbx.unit_cell((5,7,9,80,100,130))
  assert abs(adptbx.debye_waller_factor_u_iso(ucell.d_star_sq((7,8,9))/4,0.025)
             - 0.006625427) < 1.e-6
  assert abs(adptbx.debye_waller_factor_u_iso(ucell, (7, 8, 9), 0.025)
             - 0.006625427) < 1.e-6
  u = (0.1, 0.2, 0.3, -0.04, 0.05, -0.06)
  assert abs(adptbx.debye_waller_factor_u_cif(ucell, (1, 2, 3), u)
             - 0.335822877927) < 1.e-6
  #
  u_star=(0.000002, 0.000001, 0.000003, 0.0000005, 0.0000004, 0.0000006)
  assert approx_equal(
    adptbx.debye_waller_factor_u_star_gradient_coefficients(h=(5,2,-1)),
    [25, 4, 1, 20, -10, -4])
  g = adptbx.debye_waller_factor_u_star_gradients(h=(5,2,-1), u_star=u_star)
  assert approx_equal(g*1.e-3,
    [-0.49289027387879081, -0.078862443820606531, -0.019715610955151633,
     -0.39431221910303266, 0.19715610955151633, 0.078862443820606531])
  assert approx_equal(
    adptbx.debye_waller_factor_u_star_curvature_coefficients(h=(5,2,-1)),
    [625,100,25,500,-250,-100,
       16,4,80,-40,-16,
         1,20,-10,-4,
           400,-200,-80,
             100,40,
               16])
  c = adptbx.debye_waller_factor_u_star_curvatures(h=(5,2,-1), u_star=u_star)
  assert approx_equal(c*1.e-6,
    [0.24323160081641265, 0.038917056130626029, 0.0097292640326565073,
     0.1945852806531301, -0.097292640326565052, -0.038917056130626029,
     0.006226728980900164, 0.001556682245225041, 0.031133644904500817,
     -0.015566822452250408, -0.006226728980900164, 0.00038917056130626025,
     0.0077834112261252041, -0.0038917056130626021, -0.001556682245225041,
     0.15566822452250412, -0.077834112261252059, -0.031133644904500817,
     0.038917056130626029, 0.015566822452250408, 0.006226728980900164])

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
  for i in range(3):
    check_eigenvalue(u, ev[i])
  for i in range(3):
    assert abs(diag[i] - ev[i]) < 1.e-4
  if (adptbx.is_positive_definite(ev)):
    es = adptbx.eigensystem(u)
    ev = list(es.values())
    ev.sort()
    for i in range(3):
      check_eigenvalue(u, ev[i])
    for i in range(3):
      assert abs(diag[i] - ev[i]) < 1.e-4
    evec = []
    for i in range(3):
      check_eigenvector(u, es.values()[i], es.vectors(i))
      evec.extend(es.vectors(i))
    return # XXX following tests disabled for the moment
           # sometimes fail if eigenvalues are very similar but not identical
    sqrt_eval = matrix.diag(flex.sqrt(flex.double(es.values())))
    evec = matrix.sqr(evec).transpose()
    sqrt_u = evec * sqrt_eval * evec.transpose()
    u_full = matrix.sym(sym_mat3=u).elems
    assert approx_equal(u_full, (sqrt_u.transpose()*sqrt_u).elems, eps=1.e-3)
    assert approx_equal(u_full, (sqrt_u*sqrt_u.transpose()).elems, eps=1.e-3)
    assert approx_equal(u_full, (sqrt_u*sqrt_u).elems, eps=1.e-3)
    sqrt_u_plus_shifts = matrix.sym(
      sym_mat3=[x + 10*(random.random()-.5) for x in sqrt_u.as_sym_mat3()])
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
  for n_equal in range(3):
    for i_trial in range(n_trials):
      if (n_equal == 0):
        diag = [random.random() for i in range(3)]
      if (n_equal == 1):
        i = random.randrange(3)
        diag[i] = random.random()
        diag[(i+1)%3] = random.random()
        diag[(i+2)%3] = diag[(i+1)%3]
      if (n_equal == 2):
        diag = [random.random()] * 3
      exercise_eigen_core(diag)

def exercise_random_traceless_symmetry_constrained_b_cart():
  from cctbx import crystal
  cs = crystal.symmetry((10,20,30,90,90,90), "P212121")
  bc = adptbx.random_traceless_symmetry_constrained_b_cart(
    crystal_symmetry = cs, u_scale=1, u_min=0.1)
  assert approx_equal(bc[0]+bc[1]+bc[2], 0)
  assert approx_equal(bc[3],0)
  assert approx_equal(bc[4],0)
  assert approx_equal(bc[5],0)

def exercise_misc():
  import libtbx.load_env
  if (not libtbx.env.has_module("iotbx")) : return
  import iotbx.pdb
  # Pair 1: 0.5 A apart, mean displacement = 0.4 A
  # Pair 2: 1.5 A apart, mean displacement = 0.6 A
  pdb_in = iotbx.pdb.input(source_info=None, lines="""
CRYST1   10.000   11.000   12.000  70.00  80.00  90.00 P 1
HETATM    1  O  AHOH A   1       4.000   5.000   3.000  1.00 12.63           O
HETATM    2  O  BHOH A   1       4.500   5.000   3.000  1.00 12.63           O
HETATM    3  O  AHOH A   2       7.000   1.000   6.000  1.00 28.42           O
HETATM    4  O  BHOH A   2       8.500   1.000   6.000  1.00 28.42           O
END""")
  xrs = pdb_in.xray_structure_simple()
  unit_cell = xrs.unit_cell()
  sc = xrs.scatterers()
  delta12 = adptbx.intersection(
    u_1=sc[0].u_iso,
    u_2=sc[1].u_iso,
    site_1=sc[0].site,
    site_2=sc[1].site,
    unit_cell=xrs.unit_cell())
  xrs.convert_to_anisotropic()
  delta12_aniso = adptbx.intersection(
    u_1=sc[0].u_star,
    u_2=sc[1].u_star,
    site_1=sc[0].site,
    site_2=sc[1].site,
    unit_cell=xrs.unit_cell())
  # XXX on certain platforms the floating-point precision fails us
  assert approx_equal(delta12_aniso, delta12, eps=0.0000000000001)
  assert approx_equal(delta12, 0.2999, eps=0.0001)
  delta34 = adptbx.intersection(
    u_1=sc[2].u_star,
    u_2=sc[3].u_star,
    site_1=sc[2].site,
    site_2=sc[3].site,
    unit_cell=xrs.unit_cell())
  assert approx_equal(delta34, -0.300094, eps=0.000001)
  delta34b = xrs.intersection_of_scatterers(2,3)
  assert (delta34b == delta34)

def run():
  exercise_misc()
  exercise_interface()
  exercise_factor_u_star_u_iso()
  exercise_debye_waller()
  exercise_grad_u_transformations()
  exercise_eigen()
  exercise_random_traceless_symmetry_constrained_b_cart()
  print("OK")

if (__name__ == "__main__"):
  run()
