from __future__ import absolute_import, division, print_function

""" non_linear_ls_with_fixed_scale_factor, checked against known answers.

The separable class solves for the scale factor by variable projection. This
one does not: the scale is given. That is what maximum-likelihood refinement
needs, because there the scale already sits inside the likelihood's alpha and
beta, and a second one solved for on top of it does not converge to the maximum
of the likelihood.

Everything here is a linear problem, where the least-squares answer is known in
closed form, so the accumulator is compared with arithmetic rather than with
another accumulator.
"""

from scitbx.array_family import flex
from scitbx.lstbx import normal_eqns
from libtbx.test_utils import approx_equal


def design(m, n):
  """ A little design matrix and a right-hand side, nothing degenerate. """
  a = flex.double(flex.grid(m, n))
  for i in range(m):
    for j in range(n):
      a[i, j] = 1. + ((i + 1)*(j + 2) % 7)*0.37 - (i % 3)*0.11
  yo = flex.double([2.5 - 0.3*i + (i % 4)*0.9 for i in range(m)])
  w = flex.double([0.5 + (i % 5)*0.25 for i in range(m)])
  return a, yo, w


def solve_by_hand(a, yo, w, k):
  """ (k^2 A^T W A) x = k A^T W yo, formed directly. """
  m, n = a.focus()
  ata = flex.double(flex.grid(n, n), 0.)
  atb = flex.double(n, 0.)
  for i in range(m):
    for j in range(n):
      atb[j] += w[i]*yo[i]*a[i, j]*k
      for l in range(n):
        ata[j, l] += w[i]*a[i, j]*a[i, l]*k*k
  return ata, atb


def build(k, split=False):
  """ Accumulate y_c = 0 with grad y_c = the rows of A, so that the equations
  are the linear ones and the step is the whole answer. """
  m, n = 9, 3
  a, yo, w = design(m, n)
  eqns = normal_eqns.non_linear_ls_with_fixed_scale_factor(n_parameters=n)
  eqns.set_scale_factor(k)
  if not split:
    for i in range(m):
      eqns.add_equation(y_calc=0., grad_y_calc=a[i*n:(i + 1)*n],
                        y_obs=yo[i], weight=w[i])
  else:
    other = normal_eqns.non_linear_ls_with_fixed_scale_factor(n_parameters=n)
    other.set_scale_factor(k)
    for i in range(m):
      target = eqns if i < m//2 else other
      target.add_equation(y_calc=0., grad_y_calc=a[i*n:(i + 1)*n],
                          y_obs=yo[i], weight=w[i])
    eqns += other
  eqns.finalise()
  return eqns, a, yo, w, n


def exercise_normal_equations_are_the_linear_ones():
  for k in (1., 0.75, 2.5):
    eqns, a, yo, w, n = build(k)
    ata, atb = solve_by_hand(a, yo, w, k)
    b = eqns.step_equations().right_hand_side()
    assert approx_equal(list(b), list(atb), eps=1e-12), (k, list(b), list(atb))
    packed = eqns.step_equations().normal_matrix_packed_u()
    p = 0
    for i in range(n):
      for j in range(i, n):
        assert approx_equal(packed[p], ata[i, j], eps=1e-12), (k, i, j)
        p += 1
  print("\tthe normal equations are the weighted linear ones")


def exercise_objective_is_half_the_weighted_sum_of_squares():
  for k in (1., 0.75, 2.5):
    eqns, a, yo, w, n = build(k)
    # y_calc is zero, so the residual is just y_obs
    expected = 0.5*sum(w[i]*yo[i]*yo[i] for i in range(yo.size()))
    assert approx_equal(eqns.objective(), expected, eps=1e-12), (
      k, eqns.objective(), expected)
  print("\tthe objective is half the weighted sum of squares")


def exercise_the_scale_factor_is_the_one_given():
  for k in (1., 0.75, 2.5):
    eqns, a, yo, w, n = build(k)
    assert approx_equal(eqns.optimal_scale_factor(), k, eps=1e-15)
  print("\tthe scale factor is the one it was given, not one it solved for")


def exercise_merging_matches_accumulating_in_one_go():
  """ Threads accumulate into private copies and are merged with +=. """
  for k in (1., 2.5):
    whole, a, yo, w, n = build(k)
    halves, _, _, _, _ = build(k, split=True)
    assert approx_equal(list(whole.step_equations().right_hand_side()),
                        list(halves.step_equations().right_hand_side()),
                        eps=1e-12)
    assert approx_equal(list(whole.step_equations().normal_matrix_packed_u()),
                        list(halves.step_equations().normal_matrix_packed_u()),
                        eps=1e-12)
    assert approx_equal(whole.objective(), halves.objective(), eps=1e-12)
  print("\tmerging two halves matches accumulating in one go")


def exercise_solving_recovers_the_least_squares_answer():
  """ With k = 1 the step is the weighted least-squares solution itself. """
  eqns, a, yo, w, n = build(1.)
  eqns.solve()
  x = eqns.step_equations().solution()
  ata, atb = solve_by_hand(a, yo, w, 1.)
  # residual of the normal equations at the solution must vanish
  for j in range(n):
    lhs = sum(ata[j, l]*x[l] for l in range(n))
    assert approx_equal(lhs, atb[j], eps=1e-9), (j, lhs, atb[j])
  print("\tsolving recovers the weighted least-squares answer")


def run():
  exercise_normal_equations_are_the_linear_ones()
  exercise_objective_is_half_the_weighted_sum_of_squares()
  exercise_the_scale_factor_is_the_one_given()
  exercise_merging_matches_accumulating_in_one_go()
  exercise_solving_recovers_the_least_squares_answer()
  print("OK")


if __name__ == '__main__':
  run()
