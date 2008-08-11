from cctbx import adptbx
from cctbx.development import random_structure
from cctbx.development import debug_utils
from scitbx.math import full_pivoting_6_x_6
from scitbx.array_family import flex
import random
import sys
import math


def f(u_cart):
  """ A non-linear function of the 6 components u_cart """
  result = 0
  u = u_cart[0]
  for i in xrange(1,6):
    if i%2: v = u*u
    else:   v = u*u*u
    result += v * u_cart[i]
    u = u_cart[i]
  result += u*u*u * u_cart[0]
  return result


def exercise_through_space_group(flags, space_group_info):
  """ Set-up crystal and constraints """
  crystal = space_group_info.any_compatible_crystal_symmetry(
    volume=1000)
  space_group = space_group_info.group()
  unit_cell = crystal.unit_cell()
  u_star_constraints = space_group.adp_constraints()
  u_cart_constraints = space_group.cartesian_adp_constraints(unit_cell)

  """ Compatibility of all_params for u* and u_cart """
  n = u_cart_constraints.n_independent_params()
  assert n == u_star_constraints.n_independent_params()
  basis = []
  for i in xrange(n):
    v = [0]*n
    v[i] = 1
    basis.append(v)
  u_cart_basis = [ u_cart_constraints.all_params(v) for v in basis ]
  u_star_basis = [ u_star_constraints.all_params(v) for v in basis ]
  u_cart_basis_bis = [ adptbx.u_star_as_u_cart(unit_cell, u)
                       for u in u_star_basis ]
  work = flex.double(u_cart_basis)
  u_cart_basis_echelon = full_pivoting_6_x_6(work, 1e-9)
  # the vector subspaces spanned respectively by u_cart_basis and
  # by u_cart_basis_bis should be equal
  for u in u_cart_basis_bis:
    assert u_cart_basis_echelon.is_in_row_span(u, 1e-9)

  """ Test the independent gradient computation """
  eps = 1e-4
  v0 = tuple([ random.random() for i in xrange(n) ])
  u_cart_0 = u_cart_constraints.all_params(v0)
  grad_f_wrt_independent_u_cart = []
  for i in xrange(n):
    v_up = list(v0)
    v_up[i] += eps
    v_down = list(v0)
    v_down[i] -= eps
    der = (   f( u_cart_constraints.all_params(v_up)   )
            - f( u_cart_constraints.all_params(v_down) )
          ) / (2 * eps)
    grad_f_wrt_independent_u_cart.append(der)
  grad_f_wrt_u_cart = []
  for i in xrange(6):
    u_cart_up = list(u_cart_0)
    u_cart_up[i] += eps
    u_cart_down = list(u_cart_0)
    u_cart_down[i] -= eps
    der = ( f(u_cart_up) - f(u_cart_down) ) / (2 * eps)
    grad_f_wrt_u_cart.append(der)
  grad_f_wrt_independent_u_cart_1 = u_cart_constraints.independent_gradients(
    tuple(grad_f_wrt_u_cart))
  assert flex.max( flex.abs( flex.double(grad_f_wrt_independent_u_cart_1)
           - flex.double(grad_f_wrt_independent_u_cart) ) ) < 5*eps**2


def exercise_through_site_symmetry():
  from cctbx import uctbx
  uc = uctbx.unit_cell((12,12,15,90,90,120))
  from cctbx import sgtbx
  sg = sgtbx.space_group_info("P6").group()
  site_symmetry = sgtbx.site_symmetry(uc, sg, (1/3., 2/3., 0.))
  spc = site_symmetry.cartesian_adp_constraints(uc)
  sps = site_symmetry.adp_constraints()
  assert spc.n_independent_params() == sps.n_independent_params()
  site_symmetry = sgtbx.site_symmetry(uc, sg, (0.1, 0.2, 0.3))
  spc = site_symmetry.cartesian_adp_constraints(uc)
  sps = site_symmetry.adp_constraints()
  assert spc.n_independent_params() == sps.n_independent_params()


def run():
  exercise_through_site_symmetry()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:],
                                              exercise_through_space_group)
  print 'OK'

if (__name__ == "__main__"):
  run()
