from __future__ import division
from math import pi, sin, cos, asin, sqrt
import pickle
from cctbx.array_family import flex
from cctbx import uctbx
from scitbx import matrix
from libtbx.test_utils import approx_equal, not_approx_equal, show_diff
import random
import sys

def exercise_functions():
  d_star_sq_s = 1.2345
  d_star_sq_a = flex.double((1.2345, 0.1234))
  two_theta_s = 0.61725
  two_theta_a = flex.double((0.61725,0.65432))
  # forward conversions
  assert approx_equal(
    uctbx.d_star_sq_as_stol_sq(d_star_sq_s), d_star_sq_s / 4)
  assert approx_equal(
    uctbx.d_star_sq_as_stol_sq(d_star_sq_a),
    [uctbx.d_star_sq_as_stol_sq(i) for i in d_star_sq_a])
  assert approx_equal(
    uctbx.d_star_sq_as_two_stol(d_star_sq_s)**2, d_star_sq_s)
  assert approx_equal(
    uctbx.d_star_sq_as_two_stol(d_star_sq_a),
    [uctbx.d_star_sq_as_two_stol(i) for i in d_star_sq_a])
  assert approx_equal(
    uctbx.d_star_sq_as_stol(d_star_sq_s)**2, d_star_sq_s / 4)
  assert approx_equal(
    uctbx.d_star_sq_as_stol(d_star_sq_a),
    [uctbx.d_star_sq_as_stol(i) for i in d_star_sq_a])
  assert approx_equal(
    1/(uctbx.d_star_sq_as_d(d_star_sq_s)**2), d_star_sq_s)
  assert approx_equal(
    uctbx.d_star_sq_as_d(d_star_sq_a),
    [uctbx.d_star_sq_as_d(i) for i in d_star_sq_a])
  assert approx_equal(
    uctbx.d_star_sq_as_two_theta(d_star_sq_s, 1.5),
    2 * asin(1.5/2*sqrt(d_star_sq_s)))
  assert approx_equal(
    uctbx.d_star_sq_as_two_theta(d_star_sq_s, 1.5),
    uctbx.d_star_sq_as_two_theta(d_star_sq_s, 1.5, False))
  assert approx_equal(
    uctbx.d_star_sq_as_two_theta(d_star_sq_a, 1.5),
    [uctbx.d_star_sq_as_two_theta(i, 1.5) for i in d_star_sq_a])
  assert approx_equal(
    uctbx.d_star_sq_as_two_theta(d_star_sq_s, 1.5)*180/pi,
    uctbx.d_star_sq_as_two_theta(d_star_sq_s, 1.5, True))
  # reverse conversions
  for d_star_sq, two_theta in zip(
    (d_star_sq_s, d_star_sq_a),(two_theta_s, two_theta_a)):
    assert approx_equal(
      uctbx.stol_sq_as_d_star_sq(
        uctbx.d_star_sq_as_stol_sq(d_star_sq)), d_star_sq)
    assert approx_equal(
      uctbx.two_stol_as_d_star_sq(
        uctbx.d_star_sq_as_two_stol(d_star_sq)), d_star_sq)
    assert approx_equal(
      uctbx.stol_as_d_star_sq(
        uctbx.d_star_sq_as_stol(d_star_sq)), d_star_sq)
    assert approx_equal(
      uctbx.d_as_d_star_sq(
        uctbx.d_star_sq_as_d(d_star_sq)), d_star_sq)
    assert approx_equal(
      uctbx.two_theta_as_d_star_sq(
        uctbx.d_star_sq_as_two_theta(d_star_sq, 1.5), 1.5), d_star_sq)
    assert approx_equal(
      uctbx.two_theta_as_d_star_sq(two_theta, 1.5),
      uctbx.two_theta_as_d_star_sq(two_theta, 1.5, False))
    assert approx_equal(
      uctbx.two_theta_as_d_star_sq(two_theta, 1.5),
      uctbx.two_theta_as_d_star_sq(two_theta*180/pi, 1.5, True))
    assert approx_equal(
      uctbx.two_theta_as_d(two_theta, 1.5),
      uctbx.d_star_sq_as_d(uctbx.two_theta_as_d_star_sq(two_theta, 1.5)))
    assert approx_equal(
      uctbx.two_theta_as_d(two_theta, 1.5, True),
      uctbx.d_star_sq_as_d(uctbx.two_theta_as_d_star_sq(two_theta, 1.5, True)))
  #
  assert uctbx.fractional_unit_shifts(distance_frac=[0,0,0]) == (0,0,0)
  assert uctbx.fractional_unit_shifts([0.6,7.4,-0.4]) == (1,7,0)
  assert uctbx.fractional_unit_shifts([-6,3,-0.6]) == (-6,3,-1)
  site_frac_1 = [0.3,-8.6,2.1]
  for eps,expected_u2 in [(-1.e-5, 1), (1.e-5, 2)]:
    site_frac_2 = [-3,5.6,0.6-eps]
    assert uctbx.fractional_unit_shifts(
      site_frac_1=site_frac_1,
      site_frac_2=site_frac_2) == (3, -14, expected_u2)

def exercise_basic():
  d = (1,1,1,90,90,90)
  u = uctbx.unit_cell()
  assert approx_equal(u.parameters(), d)
  u = uctbx.unit_cell(d)
  assert u.parameters() == d
  assert approx_equal(u.parameters(), u.reciprocal_parameters())
  assert approx_equal(u.volume(), 1)
  assert approx_equal(u.longest_vector_sq(), 3)
  assert approx_equal(u.shortest_vector_sq(), 1)
  p = (2,3,4,80,100,110)
  for i in xrange(7):
    u = uctbx.unit_cell(p[:i])
    assert u.parameters() == p[:i] + d[i:]
    v = uctbx.unit_cell(p[:i])
    assert v.parameters() == u.parameters()
    if (i):
      assert not_approx_equal(u.parameters(), u.reciprocal_parameters())
      assert not u.is_similar_to(u.reciprocal())
      assert not u.is_similar_to(u.reciprocal(), 1.e-3)
      assert not u.is_similar_to(u.reciprocal(), 1.e-3, 1.e-3)
      assert u.is_similar_to(u.reciprocal(), 1000, 180)
    assert approx_equal(
      u.reciprocal_parameters(), u.reciprocal().parameters())
    assert approx_equal(
      u.parameters(), u.reciprocal().reciprocal_parameters())
    assert approx_equal(
      u.reciprocal_metrical_matrix(), u.reciprocal().metrical_matrix())
    assert approx_equal(
      u.metrical_matrix(), u.reciprocal().reciprocal_metrical_matrix())
    v = u.reciprocal().reciprocal()
    assert u.is_similar_to(v, 1.e-3, 1.e-3)
    assert approx_equal(u.volume(), 1/u.reciprocal().volume())
  u = uctbx.unit_cell(p)
  assert not u.is_degenerate()
  assert not u.is_degenerate(1.e-10)
  assert not u.is_degenerate(1.e-10, 1.e-5)
  assert u.is_degenerate(10)
  assert u.is_degenerate(1.e-10, 20)
  m = u.metrical_matrix()
  n = (2*2, 3*3, 4*4,
       2*3*cos(110*pi/180), 2*4*cos(100*pi/180), 3*4*cos(80*pi/180))
  assert approx_equal(m, n)
  v = uctbx.unit_cell(metrical_matrix=m)
  assert approx_equal(u.parameters(), v.parameters())
  u = uctbx.unit_cell((2,3,4))
  assert approx_equal(u.volume(), 2*3*4)
  assert approx_equal(u.longest_vector_sq(), 2*2+3*3+4*4)
  assert approx_equal(u.shortest_vector_sq(), 2*2)
  u = uctbx.unit_cell(p)
  assert approx_equal(u.volume(), 22.04006625)
  assert approx_equal(
    u.d_volume_d_params(),
    (11.020033123326023, 7.3466887488840156, 5.5100165616630115,
     0.051324088220620838, -0.051324088220620769, -0.13367230402431379))
  for alpha in xrange(70,121,10):
    for beta in xrange(70,121,10):
      for gamma in xrange(70,121,10):
        u = uctbx.unit_cell([7,11,13,alpha, beta, gamma])
        v = uctbx.unit_cell(
          orthogonalization_matrix=u.orthogonalization_matrix())
        assert v.is_similar_to(u)

def exercise_unit_cell_angles_are_feasible():
  n = 0
  n15 = 0
  for a in xrange(0,180+10,10):
    for b in xrange(0,180+10,10):
      for g in xrange(0,180+10,10):
        f = uctbx.unit_cell_angles_are_feasible(values_deg=(a,b,g))
        if (f):
          n += 1
          assert uctbx.unit_cell((10,10,10,a,b,g)).volume() > 0
        if (uctbx.unit_cell_angles_are_feasible(
                values_deg=(a,b,g), tolerance=15)):
          assert f
          n15 += 1
  assert n == 1649
  assert n15 == 1135

def exercise_frac_orth():
  u = uctbx.unit_cell(())
  assert approx_equal(
    u.fractionalization_matrix(), u.orthogonalization_matrix())
  u = uctbx.unit_cell((2,3,5))
  assert approx_equal(
    u.fractionalize((1,2,4)), (1/2., 2/3., 4/5.))
  assert approx_equal(
    u.orthogonalize((1/2., 2/3., 4/5.)), (1,2,4))
  assert approx_equal(
    u.fractionalize(flex.vec3_double([(1,2,4)])), [(1/2., 2/3., 4/5.)])
  assert approx_equal(
    u.orthogonalize(flex.vec3_double([(1/2., 2/3., 4/5.)])), [(1,2,4)])
  assert approx_equal(
    u.length((1/2.,2/3.,4/5.))**2, 1**2 + 2**2 + 4**2)
  assert approx_equal(
    u.distance((7/2.,8/3.,9/5.), (3,2,1))**2, 1**2 + 2**2 + 4**2)
  assert approx_equal(
    u.angle((0,0,0),(1/2.,0,0),(1,0,0)), 180.)
  assert approx_equal(
    u.angle((0,0,0),(1/2.,0,0),(1/2.,1/3.,0)), 90.)
  assert approx_equal(
    u.angle((1/2.,0,0),(1/2.,1/3.,0),(0,0,0)), 45.)
  assert u.angle((0,0,0),(0,0,0),(1,0,0)) is None
  assert approx_equal(
    u.dihedral((0,0,0),(1/2.,0,0),(1/2.,0,1/5.),(1,0,1/5.)), 180.)
  assert approx_equal(
    u.dihedral((0,0,0),(1/2.,0,0),(1/2.,0,1/5.),(1/2.,-1/3.,1/5.)), 90.)
  assert approx_equal(
    u.dihedral((0,0,0),(1/2.,0,0),(1/2.,0,1/5.),(1,-1/3.,1/5.)), 135.)
  assert approx_equal(
    u.dihedral((0,0,0),(1/2.,0,0),(1/2.,0,1/5.),(1,1/3.,1/5.)), -135.)
  assert u.dihedral((0,0,0),(1/2.,0,0),(1/2.,0,0),(1,0,1/5.)) is None
  assert approx_equal(
    u.mod_short_length((1/4.,2/3.,4/5.)),
    u.length((1/4.,-1/3.,-1/5.)))
  assert approx_equal(
    u.mod_short_distance((13/4.,8/3.,9/5.), (3,2,1)),
    u.length((1/4.,-1/3.,-1/5.)))
  c = flex.vec3_double(((7/2.,8/3.,9/5.), (13/4.,8/3.,9/5.)))
  assert approx_equal(
    u.min_mod_short_distance(c, (3,2,1)),
    u.mod_short_distance((13/4.,8/3.,9/5.), (3,2,1)))
  #
  u = uctbx.unit_cell((13,17,19,83,111,95))
  fm = matrix.sqr(u.fractionalization_matrix())
  assert ",".join(["%.3g" % e for e in fm.elems]) \
      == "0.0769,0.00673,0.029,0,0.059,-0.00578,0,0,0.0566"
  om = matrix.sqr(u.orthogonalization_matrix())
  assert ",".join(["%.3g" % e for e in om.elems]) \
      == "13,-1.48,-6.81,0,16.9,1.73,0,0,17.7"
  gm = matrix.sqr(u.grid_index_as_site_cart_matrix(gridding=(11,13,17)))
  pg = matrix.col((5,-7,23))
  pf = matrix.col((5/11,-7/13,23/17))
  assert approx_equal(u.orthogonalize(pf), om*pf)
  assert approx_equal(gm*pg, om*pf)
  f = flex.vec3_double(flex.random_double(size=12)*2-1)
  c = u.orthogonalize(sites_frac=f)
  assert approx_equal(u.fractionalize(sites_cart=c), f)
  for fi,ci in zip(f, c):
    assert approx_equal(u.orthogonalize(site_frac=fi), ci)
    assert approx_equal(u.fractionalize(site_cart=ci), fi)
    assert approx_equal(om*matrix.col(fi), ci)
    assert approx_equal(fm*matrix.col(ci), fi)
  #
  from cctbx import sgtbx
  s = sgtbx.rt_mx("-x,-x+y,-x+z", r_den=12)
  assert approx_equal(u.matrix_cart(rot_mx=s.r()), [
    -0.3622586, -0.1191822, -0.5137527,
    -1.435689, 0.8743934, -0.5414459,
    -1.357969, -0.1188069, 0.4878651])
  from scitbx.math import r3_rotation_axis_and_angle_from_matrix as from_matrix
  def check(u, sg):
    for s in sg:
      t = s.r().info().type()
      c = matrix.sqr(u.matrix_cart(rot_mx=s.r()))
      d = c.determinant()
      assert approx_equal(abs(d), 1)
      assert (t < 0) is (d < 0)
      fm = from_matrix(r=c*d)
      expected = {
        1: [0],
        2: [180, -180],
        3: [120, -120],
        4: [90, -90],
        6: [60, -60]}[abs(t)]
      observed = round(fm.angle(deg=True))
      if (observed not in expected):
        raise RuntimeError("%s not in %s (%s)" % (
          str(observed), str(expected), s.r().as_xyz()))
  check( # primitive settig of space group No. 230
    u=uctbx.unit_cell([10.911236359717213]*3 + [109.47122063449069]*3),
    sg=sgtbx.space_group("-I 4bd 2c 3 (y+z,x+z,x+y)"))
  check( # P 6/m m m
    u=uctbx.unit_cell((13,13,17,90,90,120)),
    sg=sgtbx.space_group("-P 6 2"))

def exercise_distance_mod_1():
  uc = uctbx.unit_cell((9,10,12,85,95,100))
  dm1 = uc.distance_mod_1(
    site_frac_1=(0.2, 0.1, -0.3),
    site_frac_2=(3.21, -6.88, 1.67))
  assert approx_equal(dm1.diff_raw, (-3.01, 6.98, -1.97))
  assert approx_equal(dm1.diff_mod, (-0.01, -0.02, 0.03))
  assert approx_equal(dm1.dist_sq, 0.1645459)
  assert dm1.unit_shifts() == (3, -7, 2)

def exercise_change_basis():
  u = uctbx.unit_cell(())
  assert approx_equal(
    u.parameters(),
    u.change_basis((1,0,0, 0,1,0, 0,0,1)).parameters())
  assert approx_equal(
    u.parameters(),
    u.change_basis((2,0,0, 0,2,0, 0,0,2), 2).parameters())
  u = uctbx.unit_cell((2,3,5))
  assert approx_equal(
    u.change_basis((0,1,0, 0,0,1, 1,0,0)).parameters(),
    (5,2,3,90,90,90))
  #
  from cctbx import sgtbx
  cb_op = sgtbx.change_of_basis_op("y,z,x").inverse()
  assert approx_equal(
    u.change_basis(cb_op=cb_op).parameters(),
    (5,2,3,90,90,90))

def exercise_miller_index_methods():
  u = uctbx.unit_cell((2,3,5))
  assert u.max_miller_indices(0.5) == (4,6,10)
  assert u.max_miller_indices(0.5, 1.e-3) == (4,6,10)
  assert approx_equal(
    u.d_star_sq((1,0,0)), 1/4.)
  assert approx_equal(
    u.d_star_sq((0,1,0)), 1/9.)
  assert approx_equal(
    u.d_star_sq((0,0,1)), 1/25.)
  u = uctbx.unit_cell((2,3,5,80,100,110))
  h = (1,2,3)
  d_star_sq_123 = 1.39498933203
  assert approx_equal(
    u.d_star_sq(h), d_star_sq_123)
  assert approx_equal(
    u.stol_sq(h), d_star_sq_123 / 4)
  assert approx_equal(
    u.two_stol(h)**2, d_star_sq_123)
  assert approx_equal(
    u.stol(h)**2, d_star_sq_123 / 4)
  assert approx_equal(
    1/(u.d(h)**2), d_star_sq_123)
  assert approx_equal(
    u.two_theta(h, 1.5), 2 * asin(1.5/2*sqrt(d_star_sq_123)))
  assert approx_equal(
    u.two_theta(h, 1.5),
    u.two_theta(h, 1.5, False))
  assert approx_equal(
    u.two_theta(h, 1.5)*180/pi,
    u.two_theta(h, 1.5, True))
  assert approx_equal(
    pow(sin(u.two_theta(h, 1.5)),2),
    u.sin_sq_two_theta(h, 1.5))
  assert approx_equal(
    sin(u.two_theta(h, 1.5)),
    u.sin_two_theta(h, 1.5))
  miller_indices = flex.miller_index(((1,2,3), (-3,4,-5), (2,-3,4)))
  for d_spacing_measure in (u.d_star_sq,
                            u.stol_sq,
                            u.two_stol,
                            u.stol,
                            u.d):
    values = d_spacing_measure(miller_indices)
    for i,v in enumerate(values):
      assert approx_equal(
        v, d_spacing_measure(miller_indices[i]))
  wavelength = 0.8
  values = u.two_theta(miller_indices, wavelength)
  for i,v in enumerate(values):
    assert approx_equal(
      v, u.two_theta(miller_indices[i], wavelength))
  for deg in (False,True):
    values = u.two_theta(miller_indices, wavelength, deg)
    for i,v in enumerate(values):
      assert approx_equal(
        v, u.two_theta(miller_indices[i], wavelength, deg))
  assert approx_equal(
    u.max_d_star_sq(miller_indices),
    u.d_star_sq((-3,4,-5)),
    eps=1e-10)
  assert approx_equal(
    u.min_max_d_star_sq(miller_indices),
    [u.d_star_sq((1,2,3)), u.d_star_sq((-3,4,-5))],
    eps=1e-10)
  values = u.sin_two_theta(miller_indices, wavelength)
  for i, v in enumerate(values):
    assert approx_equal(
      sin(u.two_theta(miller_indices[i], wavelength)),
      v)
  values = u.sin_sq_two_theta(miller_indices, wavelength)
  for i, v in enumerate(values):
    assert approx_equal(
      pow(sin(u.two_theta(miller_indices[i], wavelength)), 2),
      v)
  u = uctbx.unit_cell((2,3,5))
  rcv = u.reciprocal_space_vector((1,1,1))
  assert approx_equal(rcv, (0.5, 1/3., 0.2))
  rcvs = u.reciprocal_space_vector(miller_indices)
  assert approx_equal(rcvs,
    [(0.5, 2/3., 0.6), (-1.5, 4/3., -1.0), (1.0, -1.0, 0.8)])

def exercise_debye_waller_factor():
  u = uctbx.unit_cell((3,4,5,85,95,105))
  dw = u.debye_waller_factor
  h = (1,2,3)
  assert approx_equal(dw(miller_index=h, u_iso=0.01), 0.848878180759)
  assert approx_equal(dw(h, b_iso=1.01), 0.810924506935)
  assert approx_equal(dw(h, u_cart=[0.06,0.04,0.05,0.01,0.02,0.03]),
    0.239382185855)
  assert approx_equal(dw(h, b_cart=[1.06,1.04,1.05,1.01,1.02,1.03]),
    0.546010240906)
  assert approx_equal(dw(h, u_cif=[0.04,0.06,0.05,0.01,0.02,0.03]),
    0.251706371444)
  assert approx_equal(dw(h, u_star=[0.0004,0.0006,0.0005,0.0001,0.0002,0.0003]),
    0.781343730547)
  #
  dw(h, b_iso=-250, exp_arg_limit=60)
  dw(h, b_iso=-250, truncate_exp_arg=True)
  try: dw(h, b_iso=-250)
  except RuntimeError, e:
    assert not show_diff(str(e),
      "cctbx::adptbx::debye_waller_factor_exp:"
      " arg_limit exceeded (isotropic): arg = 51.8763 arg_limit = 50")
  else: raise Exception_expected
  try: dw(h, b_cart=[-240,-240,-240,0,0,0], exp_arg_limit=40)
  except RuntimeError, e:
    assert not show_diff(str(e),
      "cctbx::adptbx::debye_waller_factor_exp:"
      " arg_limit exceeded (anisotropic): arg = 49.8013 arg_limit = 40")
  else: raise Exception_expected
  #
  mi = flex.miller_index(((1,2,3), (-2,1,-3)))
  assert approx_equal(
    dw(miller_indices=mi, u_iso=0.01),
    [0.8488781807585344, 0.8378216424772783])
  assert approx_equal(
    dw(miller_indices=mi, u_cart=[0.06,0.04,0.05,0.01,0.02,0.03]),
    [0.2393821858545768, 0.2899416042036387])

def exercise_compare():
  u1 = uctbx.unit_cell((3,2,5,90,100,90))
  u2 = uctbx.unit_cell((2,3,5,90,80,90))
  assert u1.compare_orthorhombic(other=u1) == 0
  assert u2.compare_orthorhombic(other=u2) == 0
  assert u1.compare_orthorhombic(other=u2) == 1
  assert u2.compare_orthorhombic(other=u1) == -1
  assert u1.compare_monoclinic(
    other=u1, unique_axis=1, angular_tolerance=3) == 0
  assert u2.compare_monoclinic(
    other=u2, unique_axis=1, angular_tolerance=3) == 0
  assert u1.compare_monoclinic(
    other=u2, unique_axis=1, angular_tolerance=3) == -1
  assert u2.compare_monoclinic(
    other=u1, unique_axis=1, angular_tolerance=3) == 1
  #
  u = uctbx.unit_cell((31.8764, 6.35, 22.54, 90, 135, 90))
  c = u.change_of_basis_op_for_best_monoclinic_beta()
  assert str(c) == "a+c,b,c"
  assert c.c().r().den() == 12
  assert c.c().t().den() == 144
  u = uctbx.unit_cell((6.35, 31.8764, 16.2514, 90, 101.266, 90))
  c = u.change_of_basis_op_for_best_monoclinic_beta()
  assert str(c) == "a,b,c"
  assert c.c().r().den() == 12
  assert c.c().t().den() == 144

def exercise_is_conventional_basis():
  u = uctbx.unit_cell
  def hex(s): return u(s).is_conventional_hexagonal_basis()
  assert hex("10 10 11 90 90 120")
  assert     not hex("10 11 11 90 90 120")
  assert not hex("10 10 11 91 90 120")
  assert not hex("10 10 11 90 89 120")
  assert not hex("10 10 11 90 90 121")
  def rho(s): return u(s).is_conventional_rhombohedral_basis()
  assert     rho("10 10 10 80 80 80")
  assert not rho("10 11 10 80 80 80")
  assert not rho("10 10 11 80 80 80")
  assert not rho("10 10 10 81 80 80")
  assert not rho("10 10 10 80 81 80")
  assert not rho("10 10 10 80 80 81")

def exercise_pickle():
  u = uctbx.unit_cell((2,3,5,80,100,110))
  p = pickle.dumps(u)
  v = pickle.loads(p)
  assert u.parameters() == v.parameters()

def exercise_exceptions():
  if ("--skip" in sys.argv[1:]):
    print "SKIPPING: exercise_exceptions"
    return
  try:
    u = uctbx.unit_cell((0,0,0,0,0,0))
  except RuntimeError, e:
    assert str(e) == "cctbx Error: Unit cell parameter is zero or negative."
  else:
    raise AssertionError, 'exception expected'
  try:
    u = uctbx.unit_cell(metrical_matrix=(0,0,0,0,0,0))
  except RuntimeError, e:
    assert str(e) == "cctbx Error: Corrupt metrical matrix."
  else:
    raise AssertionError, 'exception expected'
  u = uctbx.unit_cell((2,3,5,80,100,110))
  try:
    u.two_theta((-3,4,-5), 1.5)
  except RuntimeError, e:
    assert str(e).endswith("CCTBX_ASSERT(sin_theta <= 1.0) failure.")
  else:
    raise AssertionError, 'exception expected'

def exercise_fast_minimum_reduction():
  mr = uctbx.fast_minimum_reduction(uctbx.unit_cell((1,1,1,90,90,90)))
  assert mr.iteration_limit() == 100
  assert mr.multiplier_significant_change_test() == 16
  assert mr.min_n_no_significant_change() == 2
  mr = uctbx.fast_minimum_reduction(uctbx.unit_cell((1,1,1,90,90,90)), 90)
  assert mr.iteration_limit() == 90
  assert mr.multiplier_significant_change_test() == 16
  assert mr.min_n_no_significant_change() == 2
  mr = uctbx.fast_minimum_reduction(uctbx.unit_cell((1,1,1,90,90,90)), 90,8)
  assert mr.iteration_limit() == 90
  assert mr.multiplier_significant_change_test() == 8
  assert mr.min_n_no_significant_change() == 2
  mr = uctbx.fast_minimum_reduction(uctbx.unit_cell((1,1,1,90,90,90)), 90,8,4)
  assert mr.iteration_limit() == 90
  assert mr.multiplier_significant_change_test() == 8
  assert mr.min_n_no_significant_change() == 4
  mr = uctbx.fast_minimum_reduction(uctbx.unit_cell((2,3,5,80,90,100)))
  assert approx_equal(mr.as_gruber_matrix(),(4,9,25,-5.209445,0,-2.083778))
  assert approx_equal(mr.as_niggli_matrix(),(4,9,25,-5.209445/2,0,-2.083778/2))
  assert approx_equal(mr.as_sym_mat3(),(4,9,25,-2.083778/2,0,-5.209445/2))
  assert mr.as_unit_cell().is_similar_to(uctbx.unit_cell((2,3,5,100,90,100)))
  assert approx_equal(mr.r_inv(), (-1,0,0,0,-1,0,0,0,1))
  assert mr.n_iterations() == 1
  assert not mr.termination_due_to_significant_change_test()
  assert mr.type() == 2
  mr = uctbx.fast_minimum_reduction(uctbx.unit_cell((5,3,2,50,120,130)), 8)
  assert mr.n_iterations() == 8
  assert not mr.termination_due_to_significant_change_test()
  try:
    uctbx.fast_minimum_reduction(uctbx.unit_cell((5,3,2,50,120,130)), 2, 7)
  except RuntimeError, e:
    assert str(e) == "cctbx Error: Iteration limit exceeded."
  else:
    raise AssertionError, 'exception expected'
  try:
    u = uctbx.unit_cell((2,3,5,70,120,50))
  except Exception:
    pass
  else:
    try:
      uctbx.fast_minimum_reduction(u)
    except RuntimeError, e:
      if ("--Verbose" in sys.argv[1:]):
        print "Expected:", e

class exercise_is_degenerate(object):

  def __init__(self, n_iterations=None):
    if (n_iterations is not None):
      self.n_iterations = n_iterations
    else:
      if ("--hardest" in sys.argv[1:]):
        self.n_iterations = 1000000
      elif ("--harder" in sys.argv[1:]):
        self.n_iterations = 100000
      elif ("--hard" in sys.argv[1:]):
        self.n_iterations = 10000
      else:
        self.n_iterations = 100
    self.n_stable = [0,0]
    self.n_unstable = 0
    i_iteration = 0
    rnd = random.random
    while 1:
      lengths = [rnd(), rnd(), rnd()]
      for alpha in xrange(10,180,10):
        for beta in xrange(10,180,10):
          for gamma in xrange(10,180,10):
            try:
              u = uctbx.unit_cell((2,3,5,alpha,beta,gamma))
            except Exception:
              pass
            else:
              is_degenerate = u.is_degenerate(1.e-10, 1.e-5)
              try:
                uctbx.fast_minimum_reduction(u)
                self.n_stable[int(is_degenerate)] += 1
              except RuntimeError, e:
                assert is_degenerate
                self.n_unstable += 1
              i_iteration += 1
              if (i_iteration == self.n_iterations):
                return

  def report(self):
    print "exercise_is_degenerate:"
    s = self.n_stable[0] + self.n_stable[1]
    n = self.n_iterations*0.01
    print "  n_stable:", s, self.n_stable, "= %.3g%%" % (s/n)
    print "  n_unstable:", self.n_unstable, "= %.3g%%" % (self.n_unstable/n)

def exercise_similarity_transformations():
  reference = uctbx.unit_cell(
    (79.61, 86.07, 96.9986, 89.3203, 65.7721, 62.4533))
  others = [uctbx.unit_cell(params) for params in [
    (79, 85.6519, 97.0483, 89.6713, 65.9826, 62.5374),
    (79, 85.6519, 97.0483, 68.3049, 65.9826, 62.5374)]]
  expected_transformations = [
    [(1, 0, 0, 0, 1, 0, 0, 0, 1), (1, 1, 1, 0, -1, 0, 0, 0, -1)],
    [(-1, 0, -1, 0, -1, 0, 0, 0, 1), (-1, -1, 0, 0, 1, 0, 0, 0, -1)]
  ]
  for other,expected in zip(others, expected_transformations):
    transformations = reference.similarity_transformations(
      other=other,
      relative_length_tolerance=0.02,
      absolute_angle_tolerance=2,
      unimodular_generator_range=1)
    assert list(transformations) == expected

def unit_cell_bases_mean_square_difference(self, other):
  diff_sqs = flex.double()
  for basis_vector in [(1,0,0),(0,1,0),(0,0,1)]:
    self_v = matrix.col(self.orthogonalize(basis_vector))
    other_v = matrix.col(other.orthogonalize(basis_vector))
    diff_sqs.append((self_v - other_v).norm_sq())
  return flex.mean(diff_sqs)

def exercise_bases_rmsd():
  unit_cells = [uctbx.unit_cell(params) for params in [
    (79, 85.6519, 97.0483, 89.6713, 65.9826, 62.5374),
    (79, 85.6519, 97.0483, 68.3049, 65.9826, 62.5374),
    (110, 58.4, 69.2, 90, 127, 90),
    (56.48, 56.48, 182.39, 90, 90, 120)]]
  for cell_a in unit_cells:
    for cell_b in unit_cells:
      v_cpp = cell_a.bases_mean_square_difference(cell_b)
      v_py = unit_cell_bases_mean_square_difference(cell_a, cell_b)
      assert approx_equal(v_cpp, v_py)
      if (cell_a is cell_b):
        assert approx_equal(v_cpp, 0)

def exercise_box_frac_around_sites():
  unit_cell = uctbx.unit_cell((10,10,10,90,90,120))
  buffer = 2
  sites_frac = flex.vec3_double([
    (1/2., 2/3., 0.),
    (1/2., 1/3., 0.)])
  min_, max_ = unit_cell.box_frac_around_sites(
    sites_frac=sites_frac, buffer=buffer)
  assert approx_equal(min_, (0.26905989232414967, 0.10239322565748302, -0.2))
  assert approx_equal(max_, (0.73094010767585038, 0.8976067743425169, 0.2))
  sites_cart = unit_cell.orthogonalize(sites_frac=sites_frac)
  min_, max_ = sites_cart.min(), sites_cart.max()
  min_ = unit_cell.fractionalize([m-buffer for m in min_])
  max_ = unit_cell.fractionalize([m+buffer for m in max_])
  assert approx_equal(min_, (0.017863279495408259, 0.10239322565748302, -0.2))
  assert approx_equal(max_, (0.98213672050459189, 0.8976067743425169, 0.2))
  unit_cells = [uctbx.unit_cell(params) for params in [
    (10, 15, 20, 90,  90,  90),
    (10, 10, 20, 90,  90, 120),
    (10, 10, 10, 60,  60,  60)]]
  sites_cart = flex.vec3_double([
    (2.23474, 8.72834, 4.70562),
    (3.72656, 3.28621, 9.19121),
    (-6.83519, -7.5707, 4.62386)])
  c_inv_rs = [(1,0,0, 0,1,0, 0,0,1),
              (0,1,0, 0,0,1, 1,0,0),
              (0,0,1, 1,0,0, 0,1,0)]
  for unit_cell_0 in unit_cells:
    for c_inv_r in c_inv_rs:
      unit_cell = unit_cell_0.change_basis(c_inv_r)
      sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
      min0, max0 = unit_cell.box_frac_around_sites(
        sites_cart=sites_cart)
      for x,y in zip(min0, max0): assert x < y
      for buffer in [None, 0]:
        min_, max_ = unit_cell.box_frac_around_sites(
          sites_cart=sites_cart, buffer=buffer)
        assert approx_equal(min_, min0)
        assert approx_equal(max_, max0)
        min_, max_ = unit_cell.box_frac_around_sites(
          sites_frac=sites_frac, buffer=buffer)
        assert approx_equal(min_, min0)
        assert approx_equal(max_, max0)
      for buffer in [0,3,5,7]:
        min0, max0 = unit_cell.box_frac_around_sites(
          sites_cart=sites_cart, buffer=buffer)
        for x,y in zip(min0, max0): assert x < y
        min_, max_ = unit_cell.box_frac_around_sites(
          sites_frac=sites_frac, buffer=buffer)
        assert approx_equal(min_, min0)
        assert approx_equal(max_, max0)
        min_, max_ = sites_cart.min(), sites_cart.max()
        min_ = unit_cell.fractionalize([m-buffer for m in min_])
        max_ = unit_cell.fractionalize([m+buffer for m in max_])
        if (unit_cell_0 is unit_cells[0]):
          assert approx_equal(min_, min0)
          assert approx_equal(max_, max0)
        elif (buffer == 0 and unit_cell_0 is unit_cells[1]):
          assert approx_equal(min_, min0)
          if (c_inv_r is c_inv_rs[2]):
            assert approx_equal(max_, max0)
          else:
            assert not_approx_equal(max_, max0)
        else:
          assert not_approx_equal(min_, min0)
          assert not_approx_equal(max_, max0)

def exercise_non_crystallographic_unit_cell_with_the_sites_in_its_center():
  sites_cart = flex.vec3_double([(-5.,-5.,-5.)])
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
                    sites_cart   = sites_cart,
                    buffer_layer = 5)
  assert approx_equal(box.unit_cell.parameters(), (10, 10, 10, 90, 90, 90))
  assert approx_equal(box.sites_cart, [(5.0, 5.0, 5.0)])
  assert box.crystal_symmetry().space_group_info().type().number() == 1

def exercise_tensor_rank_2_orth_and_frac_linear_maps():
  from cctbx import adptbx, sgtbx
  p1 = sgtbx.space_group_info('P1')
  for i in xrange(100):
    uc = p1.any_compatible_unit_cell(27)
    u_star = matrix.col.random(n=6, a=0, b=1)
    u_iso_ref = adptbx.u_star_as_u_iso(uc, u_star)
    u_iso = matrix.col(uc.u_star_to_u_iso_linear_form()).dot(u_star)
    assert approx_equal(u_iso, u_iso_ref, eps=1e-15)
    u_cart_ref = adptbx.u_star_as_u_cart(uc, u_star)
    u_cart = matrix.sqr(uc.u_star_to_u_cart_linear_map()) * u_star
    assert approx_equal(u_cart, u_cart_ref, eps=1e-15)
    u_cif_ref = adptbx.u_star_as_u_cif(uc, u_star)
    u_cif = matrix.diag(uc.u_star_to_u_cif_linear_map())*(u_star)
    assert approx_equal(u_cif, u_cif_ref)

def exercise_d_metrical_matrix_d_params():
  def finite_differences(unit_cell, eps=1e-6):
    grads = []
    for i in range(6):
      params = list(unit_cell.parameters())
      params[i] += eps
      uc = uctbx.unit_cell(parameters=params)
      qm = matrix.col(uc.metrical_matrix())
      params[i] -= 2*eps
      uc = uctbx.unit_cell(parameters=params)
      qp = matrix.col(uc.metrical_matrix())
      dq = (qm-qp)/(2*eps)
      grads.extend(list(dq))
    grads = flex.double(grads)
    grads.resize(flex.grid((6,6)))
    return grads.matrix_transpose()
  from cctbx import sgtbx
  p1 = sgtbx.space_group_info('P1')
  uc = p1.any_compatible_unit_cell(27)
  grads = uc.d_metrical_matrix_d_params()
  fd_grads = finite_differences(uc)
  assert approx_equal(grads, fd_grads)
  sgi = sgtbx.space_group_info('I-4')
  uc = sgi.any_compatible_unit_cell(volume=18000)
  grads = uc.d_metrical_matrix_d_params()
  fd_grads = finite_differences(uc)
  assert approx_equal(grads, fd_grads)

def exercise_downstream_methods():
  uc = uctbx.unit_cell((10,10,10,90,90,90))
  assert str(uc.lattice_symmetry_group().info()) == "P 4 3 2"
  for anomalous_flag in [False, True]:
    ms = uc.complete_miller_set_with_lattice_symmetry(
      d_min=2-1e-6, anomalous_flag=anomalous_flag)
    assert str(ms.space_group_info()) == "P 4 3 2"
    if (anomalous_flag):
      assert ms.indices().size() == 28
    else:
      assert ms.indices().size() == 26

def run():
  exercise_d_metrical_matrix_d_params()
  exercise_tensor_rank_2_orth_and_frac_linear_maps()
  exercise_non_crystallographic_unit_cell_with_the_sites_in_its_center()
  exercise_functions()
  exercise_basic()
  exercise_unit_cell_angles_are_feasible()
  exercise_frac_orth()
  exercise_distance_mod_1()
  exercise_change_basis()
  exercise_miller_index_methods()
  exercise_debye_waller_factor()
  exercise_compare()
  exercise_is_conventional_basis()
  exercise_pickle()
  exercise_exceptions()
  exercise_fast_minimum_reduction()
  exercise_similarity_transformations()
  exercise_bases_rmsd()
  exercise_box_frac_around_sites()
  exercise_downstream_methods()
  e = exercise_is_degenerate()
  if (e.n_iterations > 100):
    e.report()
  print "OK"

if (__name__ == "__main__"):
  run()
