from __future__ import absolute_import, division, print_function
from cctbx import adptbx, sgtbx
from cctbx.development import debug_utils
from scitbx.math import row_echelon_full_pivoting
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
import random
import sys
from six.moves import range


def f(u_cart):
  """ A non-linear function of the 6 components u_cart """
  result = 0
  u = u_cart[0]
  for i in range(1,6):
    if i%2: v = u*u
    else:   v = u*u*u
    result += v * u_cart[i]
    u = u_cart[i]
  result += u*u*u * u_cart[0]
  return result


def exercise_all_wyckoff(flags, space_group_info):
  """ Set-up crystal and constraints """
  crystal = space_group_info.any_compatible_crystal_symmetry(
    volume=1000)
  unit_cell = crystal.unit_cell()
  wyckoffs = space_group_info.wyckoff_table()
  for i in range(wyckoffs.size()):
    wyckoff_pos = wyckoffs.position(i)
    print("%s," % wyckoff_pos.letter(), end=' ')
    special_op = wyckoff_pos.special_op()
    exact_site = eval("(%s)" % special_op, {'x':0.1, 'y':0.2, 'z':0.3})
    site_symmetry = sgtbx.site_symmetry(crystal.unit_cell(),
                                        crystal.space_group(),
                                        exact_site)
    u_star_constraints = site_symmetry.adp_constraints()
    u_cart_constraints = site_symmetry.cartesian_adp_constraints(unit_cell)

    """ Compatibility of all_params for u* and u_cart """
    n = u_cart_constraints.n_independent_params()
    assert n == u_star_constraints.n_independent_params()
    basis = []
    for i in range(n):
      v = [0]*n
      v[i] = 1
      basis.append(v)
    u_cart_basis = [ u_cart_constraints.all_params(v) for v in basis ]
    u_star_basis = [ u_star_constraints.all_params(v) for v in basis ]
    u_cart_basis_bis = [ adptbx.u_star_as_u_cart(unit_cell, u)
                         for u in u_star_basis ]
    a_work = flex.double(u_cart_basis)
    u_cart_basis_echelon = row_echelon_full_pivoting(
      a_work=a_work, min_abs_pivot=1e-9)
    # the vector subspaces spanned respectively by u_cart_basis and
    # by u_cart_basis_bis should be equal
    for u in u_cart_basis_bis:
      assert u_cart_basis_echelon.is_in_row_space(
        x=flex.double(u), epsilon=1e-9)

    """ Test the independent gradient computation """
    eps = 1e-4
    v0 = tuple([ random.random() for i in range(n) ])
    u_cart_0 = u_cart_constraints.all_params(v0)
    grad_f_wrt_independent_u_cart = []
    for i in range(n):
      v_up = list(v0)
      v_up[i] += eps
      v_down = list(v0)
      v_down[i] -= eps
      der = (   f( u_cart_constraints.all_params(v_up)   )
              - f( u_cart_constraints.all_params(v_down) )
            ) / (2 * eps)
      grad_f_wrt_independent_u_cart.append(der)
    grad_f_wrt_u_cart = []
    for i in range(6):
      u_cart_up = list(u_cart_0)
      u_cart_up[i] += eps
      u_cart_down = list(u_cart_0)
      u_cart_down[i] -= eps
      der = ( f(u_cart_up) - f(u_cart_down) ) / (2 * eps)
      grad_f_wrt_u_cart.append(der)
    grad_f_wrt_independent_u_cart_1 = (
      u_cart_constraints.independent_gradients(tuple(grad_f_wrt_u_cart)))
    assert flex.max( flex.abs( flex.double(grad_f_wrt_independent_u_cart_1)
             - flex.double(grad_f_wrt_independent_u_cart) ) ) < 5*eps**2

    """ Check independent_params """
    v = tuple([ random.random() for i in range(n) ])
    u_cart = u_cart_constraints.all_params(v)
    w = u_cart_constraints.independent_params(u_cart)
    assert approx_equal(v, w, eps=1e-12)

  print()

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:],
                                              exercise_all_wyckoff)

if (__name__ == "__main__"):
  run()
