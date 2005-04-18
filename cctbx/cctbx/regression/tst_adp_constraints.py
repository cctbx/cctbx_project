from cctbx import miller
from cctbx import crystal
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
  return -2*math.pi**2 * matrix.row([
    h[0] * h[0],
    h[1] * h[1],
    h[2] * h[2],
    2 * h[0] * h[1],
    2 * h[0] * h[2],
    2 * h[1] * h[2]]) * adptbx.debye_waller_factor_u_star(h, u_star)

def dw_analytical_grads_full_sym(space_group, h, u_star):
  result = matrix.row([0]*6)
  n = 0
  for eq in miller.sym_equiv_indices(space_group, h).p1_listing(False):
    result += matrix.row(dw_analytical_grads_full(eq.hr(), u_star))
    n += 1
  return tuple(result / n)

def dw_finite_difference_grads_full_sym(space_group, h, u_star, eps):
  result = []
  for i,u in enumerate(u_star):
    u_star_shifted = list(u_star)
    u_star_shifted[i] = u + eps
    dw_plus = dw_sym(space_group, h, u_star_shifted, no_assert=True)
    u_star_shifted[i] = u - eps
    dw_minus = dw_sym(space_group, h, u_star_shifted, no_assert=True)
    result.append((dw_plus-dw_minus) / (2*eps))
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
    u_full = list(adp_constraints.all_params(flex.double(u_shifted)))
    dw_plus = dw_sym(space_group, h, u_full, no_assert=True)
    u_shifted[i] = u - eps
    u_full = list(adp_constraints.all_params(flex.double(u_shifted)))
    dw_minus = dw_sym(space_group, h, u_full, no_assert=True)
    result.append((dw_plus-dw_minus) / (2*eps))
  return result

def show_scaled(values, scale, fmt="%10.3f"):
  for value in values:
    print fmt % (value*scale),
  print

def are_similar(a1, a2):
  a1 = flex.double(a1)
  a2 = flex.double(a2)
  m = max(flex.max(flex.abs(a1)), flex.max(flex.abs(a2)))
  if (m > 0):
    a1 /= m
    a2 /= m
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
    print list(adp_constraints.independent_flags)
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
      print "g_full_fin:  ",
      show_scaled(values=g_full_fin, scale=u_scale)
      print "g_full_ana:    ",
      show_scaled(values=g_full_ana, scale=u_scale)
    assert are_similar(g_full_fin, g_full_ana)
    g_full_ana_ave = adp_constraints.sym_gradients(asu_gradients=g0)
    if (verbose):
      print "g_full_ana_ave:",
      show_scaled(values=g_full_ana_ave, scale=u_scale)
    assert are_similar(g_full_ana_ave, g_full_ana)
    g_fin_indep = dw_finite_difference_grads_indep(
      space_group=space_group,
      adp_constraints=adp_constraints,
      h=h,
      u_star=u_star,
      eps=eps)
    g_ana_indep = adp_constraints.independent_gradients(all_gradients=g0)
    if (verbose):
      print "g_fin_indep:",
      show_scaled(values=g_fin_indep, scale=u_scale)
      print "g_ana_indep:",
      show_scaled(values=g_ana_indep, scale=u_scale)
    assert are_similar(g_ana_indep, g_fin_indep)

def run_call_back(flags, space_group_info):
  exercise(space_group_info, verbose=flags.Verbose)
  if (space_group_info.group().n_ltr() != 1):
    exercise(space_group_info.primitive_setting(), verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run()
