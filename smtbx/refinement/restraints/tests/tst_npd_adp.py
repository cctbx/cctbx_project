from __future__ import absolute_import, division, print_function

from cctbx import adp_restraints, adptbx, xray, uctbx
from cctbx.adp_restraints import adp_restraint_params
from cctbx.array_family import flex
from cctbx.xray import parameter_map
from scitbx import matrix
from smtbx.refinement import restraints
from libtbx.test_utils import approx_equal
from six.moves import range

adp = adp_restraints


def one_atom_structure(u_cart):
  """ A single anisotropic C in a deliberately triclinic cell, so the
  u_star <-> u_cart map is not the identity and the off-diagonal handling in
  linearise is actually exercised. """
  from cctbx import crystal
  uc = uctbx.unit_cell((7.1, 8.3, 9.7, 85.0, 95.0, 105.0))
  cs = crystal.symmetry(unit_cell=uc, space_group_symbol="P1")
  sc = xray.scatterer(label="C1", site=(0.11, 0.23, 0.37))
  sc.flags.set_use_u_aniso(True)
  sc.u_star = adptbx.u_cart_as_u_star(uc, u_cart)
  sc.flags.set_grad_u_aniso(True)
  xs = xray.structure(crystal_symmetry=cs)
  xs.add_scatterer(sc)
  return xs


def s_of(u_cart):
  a, b, c, d, e, f = u_cart
  u_eq = (a + b + c) / 3
  det = a*b*c + 2*d*e*f - a*f*f - b*e*e - c*d*d
  return det / u_eq**3


def exercise_s_and_activation():
  """ s is the normalised determinant, 1 for a sphere and 0 at the NPD
  boundary, and the restraint only bites below the target. """
  # a sphere
  r = adp.npd_adp((0.02, 0.02, 0.02, 0, 0, 0), weight=1e4, s_target=0.5)
  assert approx_equal(r.s(), 1.0)
  assert not r.active() and r.delta() == 0 and r.residual() == 0

  # a healthy, mildly anisotropic atom: PD, s well above a small target
  u = (0.05, 0.05, 0.005, 0.005, 0.003, 0.002)
  r = adp.npd_adp(u, weight=1e4, s_target=0.001)
  assert r.s() > 0.2
  assert not r.active()

  # the same atom against a target above its s: now restrained
  r = adp.npd_adp(u, weight=1e4, s_target=0.5)
  assert r.active()
  assert approx_equal(r.s(), s_of(u))
  assert approx_equal(r.delta(), 0.5 - s_of(u))
  assert approx_equal(r.residual(), 1e4 * (0.5 - s_of(u))**2)

  # a non-positive-definite atom (one negative eigenvalue): s < 0, restrained
  npd = (0.05, 0.05, -0.01, 0.0, 0.0, 0.0)
  r = adp.npd_adp(npd, weight=1e4, s_target=0.001)
  assert r.s() < 0 and r.active()
  print("\ts and activation OK")


def exercise_u_cart_gradient():
  """ The u_cart gradient built into the restraint (cart_grad_, hence the
  determinant math) against a finite difference of delta in u_cart. Checked
  through gradients() = d(residual)/d(u_cart) = 2 w delta d(delta)/d(u_cart). """
  u = list((0.05, 0.048, 0.006, 0.004, 0.003, 0.002))
  w, s_target = 1e4, 0.5
  r = adp.npd_adp(tuple(u), weight=w, s_target=s_target)
  assert r.active()
  g = r.gradients()  # d(residual)/d(u_cart), off-diagonals NOT doubled here
  eps = 1e-9
  for j in range(6):
    up = u[:]; up[j] += eps
    dn = u[:]; dn[j] -= eps
    # residual = w delta^2; central difference in u_cart component j
    d_res = (adp.npd_adp(tuple(up), weight=w, s_target=s_target).residual()
             - adp.npd_adp(tuple(dn), weight=w, s_target=s_target).residual()
             ) / (2*eps)
    assert approx_equal(g[j], d_res, eps=1e-3), (j, g[j], d_res)
  print("\tu_cart gradient matches finite differences")


def exercise_linearise_design_row():
  """ The least-squares path: the design-matrix row against a finite difference
  of delta in the refined u_star parameters, with the off-diagonal factor of
  two the cctbx convention applies. This is the row Olex2's CGLS consumes. """
  u_cart = (0.05, 0.048, 0.006, 0.004, 0.003, 0.002)
  w, s_target = 1e4, 0.5
  xs = one_atom_structure(u_cart)
  uc = xs.unit_cell()
  pm = xs.parameter_map()

  proxies = adp.shared_npd_adp_proxy([
    adp.npd_adp_proxy((0,), weight=w, s_target=s_target)])
  mgr = restraints.manager(npd_adp_proxies=proxies)
  eqns = mgr.build_linearised_eqns(xs, pm)
  design = eqns.design_matrix.as_dense_matrix()
  assert eqns.n_restraints() == 1
  delta0 = eqns.deltas[0]
  assert approx_equal(delta0, s_target - s_of(u_cart))

  # finite difference of delta w.r.t. each u_star component
  u_star = list(xs.scatterers().extract_u_star()[0])
  eps = 1e-9
  u_aniso = pm[0].u_aniso
  for j in range(6):
    hp = list(u_star); hp[j] += eps
    hm = list(u_star); hm[j] -= eps
    dp = s_target - s_of(adptbx.u_star_as_u_cart(uc, hp))
    dm = s_target - s_of(adptbx.u_star_as_u_cart(uc, hm))
    d_delta = (dp - dm) / (2*eps)
    if j > 2:
      d_delta *= 2   # off diagonals count twice, as tst_restraints.py does
    assert approx_equal(design[(0, u_aniso + j)], d_delta, eps=1e-2), \
      (j, design[(0, u_aniso + j)], d_delta)
  print("\tlinearise design row matches finite differences (with the x2)")


def exercise_inactive_atom_is_a_zero_row():
  """ A healthy atom still emits exactly one row -- the manager counted it --
  but with zero weight, so it adds nothing to the normal equations. """
  healthy = (0.02, 0.021, 0.019, 0.001, 0.0, 0.0)   # PD, s near 1
  xs = one_atom_structure(healthy)
  proxies = adp.shared_npd_adp_proxy([
    adp.npd_adp_proxy((0,), weight=1e4, s_target=0.001)])
  mgr = restraints.manager(npd_adp_proxies=proxies)
  eqns = mgr.build_linearised_eqns(xs, xs.parameter_map())
  assert eqns.n_restraints() == 1
  assert eqns.weights[0] == 0
  assert eqns.deltas[0] == 0
  print("\tinactive atom is a single zero-weight row")


def run():
  exercise_s_and_activation()
  exercise_u_cart_gradient()
  exercise_linearise_design_row()
  exercise_inactive_atom_is_a_zero_row()
  print("OK")


if __name__ == '__main__':
  run()
