from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import adptbx
from cctbx.development import debug_utils
from cctbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal
import scitbx.math
import math
import sys

def dw_sym(space_group, h, u_star, no_assert=False):
  result = 0
  n = 0
  for eq in miller.sym_equiv_indices(space_group, h).p1_listing(False):
    result += adptbx.debye_waller_factor_u_star(eq.hr(), u_star)
    n += 1
  if (not no_assert):
    assert approx_equal(result/n, adptbx.debye_waller_factor_u_star(h, u_star))
  return result / n

def dw_analytical_grads_full(h, u_star):
  return list(adptbx.debye_waller_factor_u_star_gradients(h, u_star))

def dw_analytical_grads_full_sym(space_group, h, u_star):
  result = matrix.row([0]*6)
  n = 0
  for eq in miller.sym_equiv_indices(space_group, h).p1_listing(False):
    result += matrix.row(dw_analytical_grads_full(eq.hr(), u_star))
    n += 1
  return tuple(result / n)

def dw_finite_difference_grads_full_sym(space_group, h, u_star, eps):
  result = []
  for i in xrange(6):
    dws = []
    for signed_eps in [eps,-eps]:
      u_star_eps = list(u_star)
      u_star_eps[i] += signed_eps
      dws.append(dw_sym(space_group, h, u_star_eps, no_assert=True))
    result.append((dws[0]-dws[1]) / (2*eps))
  return result

def dw_finite_difference_grads_indep(
      space_group,
      adp_constraints,
      h,
      u_star,
      eps):
  u_indep = adp_constraints.independent_params(u_star)
  result = []
  for i,u in enumerate(u_indep):
    u_shifted = list(u_indep)
    u_shifted[i] = u + eps
    u_full = adp_constraints.all_params(u_shifted)
    dw_plus = dw_sym(space_group, h, u_full, no_assert=True)
    u_shifted[i] = u - eps
    u_full = adp_constraints.all_params(u_shifted)
    dw_minus = dw_sym(space_group, h, u_full, no_assert=True)
    result.append((dw_plus-dw_minus) / (2*eps))
  return result

def d_dw_d_uij_finite(h, u_star, ij, eps=1.e-6):
  dws = []
  for signed_eps in [eps,-eps]:
    u_star_eps = list(u_star)
    u_star_eps[ij] += signed_eps
    dws.append(adptbx.debye_waller_factor_u_star(h, u_star_eps))
  return (dws[0]-dws[1]) / (2*eps)

def d_dw2_d_uij_d_uij_finite(h, u_star, ij1, ij2, eps=1.e-6):
  dws = []
  for signed_eps in [eps,-eps]:
    u_star_eps = list(u_star)
    u_star_eps[ij2] += signed_eps
    dws.append(d_dw_d_uij_finite(h, u_star_eps, ij=ij1))
  return (dws[0]-dws[1]) / (2*eps)

def dw_finite_difference_curvatures(h, u_star):
  result = flex.double()
  for ij1 in xrange(6):
    for ij2 in xrange(ij1,6):
      result.append(d_dw2_d_uij_d_uij_finite(h, u_star, ij1, ij2))
  return result

def d_dw_d_u_indep_finite(adp_constraints, h, u_indep, eps=1.e-6):
  result = []
  for i_indep in xrange(len(u_indep)):
    dws = []
    for signed_eps in [eps,-eps]:
      u_indep_eps = list(u_indep)
      u_indep_eps[i_indep] += signed_eps
      u_eps = adp_constraints.all_params(u_indep_eps)
      dws.append(adptbx.debye_waller_factor_u_star(h, u_eps))
    result.append((dws[0]-dws[1]) / (2*eps))
  return result

def d_dw2_d_u_indep_d_u_indep_finite(adp_constraints, h, u_indep, eps=1.e-6):
  result = []
  for i_indep in xrange(len(u_indep)):
    dws = []
    for signed_eps in [eps,-eps]:
      u_indep_eps = list(u_indep)
      u_indep_eps[i_indep] += signed_eps
      dws.append(d_dw_d_u_indep_finite(adp_constraints, h, u_indep_eps))
    row = [(x-y)/(2*eps) for x,y in zip(dws[0], dws[1])]
    result.append(row)
  return result

def dw_indep_finite_difference_curvatures(adp_constraints, h, u_star):
  u_indep = adp_constraints.independent_params(u_star)
  return flex.double(d_dw2_d_u_indep_d_u_indep_finite(
    adp_constraints, h, u_indep)).matrix_upper_diagonal()

def p2_curv(h, u_star):
  """\
h = {h0,h1,h2}
u = {{u00,0,u02},{0,u11,0},{u02,0,u22}}
dw=Exp[mtps*h.u.h]
FortranForm[D[dw,u00,u00]/dw]
FortranForm[D[dw,u00,u11]/dw]
FortranForm[D[dw,u00,u22]/dw]
FortranForm[D[dw,u00,u02]/dw]
FortranForm[D[dw,u11,u00]/dw]
FortranForm[D[dw,u11,u11]/dw]
FortranForm[D[dw,u11,u22]/dw]
FortranForm[D[dw,u11,u02]/dw]
FortranForm[D[dw,u22,u00]/dw]
FortranForm[D[dw,u22,u11]/dw]
FortranForm[D[dw,u22,u22]/dw]
FortranForm[D[dw,u22,u02]/dw]
FortranForm[D[dw,u02,u00]/dw]
FortranForm[D[dw,u02,u11]/dw]
FortranForm[D[dw,u02,u22]/dw]
FortranForm[D[dw,u02,u02]/dw]
"""
  dw = adptbx.debye_waller_factor_u_star(h, u_star)
  h0, h1, h2 = h
  u00 = u_star[0]
  u11 = u_star[1]
  u22 = u_star[2]
  u02 = u_star[4]
  mtps = -2*math.pi**2
  return [
    [dw * (h0**4*mtps**2),
     dw * (h0**2*h1**2*mtps**2),
     dw * (h0**2*h2**2*mtps**2),
     dw * (2*h0**3*h2*mtps**2)],
    [dw * (h0**2*h1**2*mtps**2),
     dw * (h1**4*mtps**2),
     dw * (h1**2*h2**2*mtps**2),
     dw * (2*h0*h1**2*h2*mtps**2)],
    [dw * (h0**2*h2**2*mtps**2),
     dw * (h1**2*h2**2*mtps**2),
     dw * (h2**4*mtps**2),
     dw * (2*h0*h2**3*mtps**2)],
    [dw * (2*h0**3*h2*mtps**2),
     dw * (2*h0*h1**2*h2*mtps**2),
     dw * (2*h0*h2**3*mtps**2),
     dw * (4*h0**2*h2**2*mtps**2)]]

def p4_curv(h, u_star):
  """\
h = {h0,h1,h2}
u = {{u11,0,0},{0,u11,0},{0,0,u22}}
dw=Exp[mtps*h.u.h]
FortranForm[D[dw,u11,u11]/dw]
FortranForm[D[dw,u11,u22]/dw]
FortranForm[D[dw,u22,u11]/dw]
FortranForm[D[dw,u22,u22]/dw]
"""
  dw = adptbx.debye_waller_factor_u_star(h, u_star)
  h0, h1, h2 = h
  u11 = u_star[1]
  u22 = u_star[2]
  mtps = -2*math.pi**2
  return [
    [dw * ((h0**2 + h1**2)**2*mtps**2),
     dw * ((h0**2 + h1**2)*h2**2*mtps**2)],
    [dw * ((h0**2 + h1**2)*h2**2*mtps**2),
     dw * (h2**4*mtps**2)]]

def p3_curv(h, u_star):
  """\
h = {h0,h1,h2}
u = {{2*u01,u01,0},{u01,2*u01,0},{0,0,u22}}
dw=Exp[mtps*h.u.h]
FortranForm[D[dw,u22,u22]/dw]
FortranForm[D[dw,u22,u01]/dw]
FortranForm[D[dw,u01,u22]/dw]
FortranForm[D[dw,u01,u01]/dw]
"""
  dw = adptbx.debye_waller_factor_u_star(h, u_star)
  h0, h1, h2 = h
  u22 = u_star[2]
  u01 = u_star[3]
  mtps = -2*math.pi**2
  return [
    [dw * (h2**4*mtps**2),
     dw * ((h0*(2*h0 + h1) + h1*(h0 + 2*h1))*h2**2*mtps**2)],
    [dw * ((h0*(2*h0 + h1) + h1*(h0 + 2*h1))*h2**2*mtps**2),
     dw * ((h0*(2*h0 + h1) + h1*(h0 + 2*h1))**2*mtps**2)]]

def p23_curv(h, u_star):
  """\
h = {h0,h1,h2}
u = {{u22,0,0},{0,u22,0},{0,0,u22}}
dw=Exp[mtps*h.u.h]
FortranForm[D[dw,u22,u22]/dw]
"""
  dw = adptbx.debye_waller_factor_u_star(h, u_star)
  h0, h1, h2 = h
  u22 = u_star[2]
  mtps = -2*math.pi**2
  return [
    [dw * ((h0**2 + h1**2 + h2**2)**2*mtps**2)]]

def show_scaled(values, scale, fmt="%10.3f"):
  for value in values:
    print fmt % (value*scale),
  print

def are_similar(a1, a2):
  a1 = flex.double(a1)
  a2 = flex.double(a2)
  m = max(flex.max(flex.abs(a1)), flex.max(flex.abs(a2)))
  if (m > 0):
    a1 = a1 / m
    a2 = a2 / m
  return approx_equal(a1, a2)

def exercise(space_group_info, verbose):
  crystal_symmetry = crystal.symmetry(
    unit_cell=space_group_info.any_compatible_unit_cell(volume=1000),
      space_group_info=space_group_info)
  space_group = space_group_info.group()
  adp_constraints = space_group.adp_constraints(
    initialize_gradient_handling=True)
  m = adp_constraints.row_echelon_form
  if (verbose):
    print matrix.rec(m, (m.size()/6, 6)).mathematica_form(
      one_row_per_line=True)
    print list(adp_constraints.independent_indices)
  u_cart_p1 = adptbx.random_u_cart()
  u_star_p1 = adptbx.u_cart_as_u_star(crystal_symmetry.unit_cell(), u_cart_p1)
  u_star = space_group.average_u_star(u_star_p1)
  u_scale = crystal_symmetry.unit_cell().volume()**(2./3.)
  eps = 1.e-10 * u_scale
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry, d_min=3, anomalous_flag=False)
  for h in miller_set.indices():
    if (verbose):
      print "dw: %8.6f" % dw_sym(space_group, h, u_star)
    g0 = dw_analytical_grads_full(h, u_star)
    ghrs = []
    for i in xrange(space_group.n_smx()):
      r = space_group(i).r()
      assert r.den() == 1
      hr = matrix.row(h) * matrix.sqr(r.num())
      ghr = dw_analytical_grads_full(hr, u_star)
      gtr = scitbx.math.tensor_rank_2_gradient_transform(a=r.num(), g=g0)
      if (verbose):
        print "ghr:",
        show_scaled(values=ghr, scale=u_scale)
        print "gtr:",
        show_scaled(values=gtr, scale=u_scale)
      assert are_similar(gtr, ghr)
      ghrs.append(ghr)
    g_full_fin = dw_finite_difference_grads_full_sym(
      space_group=space_group,
      h=h,
      u_star=u_star,
      eps=eps)
    g_full_ana = dw_analytical_grads_full_sym(
      space_group=space_group,
      h=h,
      u_star=u_star)
    if (verbose):
      print "g_full_fin:",
      show_scaled(values=g_full_fin, scale=u_scale)
      print "g_full_ana:",
      show_scaled(values=g_full_ana, scale=u_scale)
    assert are_similar(g_full_fin, g_full_ana)
    g_full_ave = adp_constraints.sym_gradients(asu_gradients=g0)
    if (verbose):
      print "g_full_ave:",
      show_scaled(values=g_full_ave, scale=u_scale)
    assert are_similar(g_full_ave, g_full_ana)
    g_fin_indep = dw_finite_difference_grads_indep(
      space_group=space_group,
      adp_constraints=adp_constraints,
      h=h,
      u_star=u_star,
      eps=eps)
    g_ave_indep = adp_constraints.independent_gradients(
      all_gradients=g_full_ave)
    if (verbose):
      print "g_fin_indep:",
      show_scaled(values=g_fin_indep, scale=u_scale)
      print "g_ave_indep:",
      show_scaled(values=g_ave_indep, scale=u_scale)
    assert are_similar(g_ave_indep, g_fin_indep)
    for ghr in ghrs:
      g_ana_indep = adp_constraints.independent_gradients(
        all_gradients=ghr)
      if (verbose):
        print "g_ana_indep:",
        show_scaled(values=g_ana_indep, scale=u_scale)
      assert are_similar(g_ana_indep, g_ave_indep)
    #
    f2 = dw_finite_difference_curvatures(h=h, u_star=u_star)
    a2 = adptbx.debye_waller_factor_u_star_curvatures(h=h, u_star=u_star)
    assert are_similar(a2, f2)
    if2 = dw_indep_finite_difference_curvatures(adp_constraints, h, u_star)
    ia2 = adp_constraints.independent_curvatures(all_curvatures=a2)
    assert are_similar(ia2, if2)
    #
    ma2 = None
    if (str(space_group_info) == "P 1 2 1"):
      assert list(adp_constraints.independent_indices) == [0,1,2,4]
      ma2 = p2_curv(h, u_star)
    elif (str(space_group_info) == "P 4"):
      assert list(adp_constraints.independent_indices) == [1,2]
      ma2 = p4_curv(h, u_star)
    elif (str(space_group_info) in ["P 3", "P 6"]):
      assert list(adp_constraints.independent_indices) == [2,3]
      ma2 = p3_curv(h, u_star)
    elif (str(space_group_info) == "P 2 3"):
      assert list(adp_constraints.independent_indices) == [2]
      ma2 = p23_curv(h, u_star)
    if (ma2 is not None):
      assert are_similar(ia2, flex.double(ma2).matrix_upper_diagonal())

def run_call_back(flags, space_group_info):
  exercise(space_group_info, verbose=flags.Verbose)
  if (space_group_info.group().n_ltr() != 1):
    exercise(space_group_info.primitive_setting(), verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run()
