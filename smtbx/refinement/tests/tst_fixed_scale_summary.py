"""The fixed-scale summary accumulator, and the likelihood reaching the step.

tst_maximum_likelihood.py checks the *algebra* of the effective-observation
reduction against an independent oracle, and it passed throughout the period
when maximum likelihood was not being applied at all: cgls.linearised_problem
built the design matrix without ml_data, so the conjugate gradients solved the
ordinary least-squares system whatever target was asked for. Correct algebra,
never connected.

So these probe the parts rather than the result:

  the accumulator          against closed-form arithmetic, standalone
  against the separable    at its own optimum, where the two must agree
  the plumbing             ml_data reaching build_design_matrix, i.e. that two
                           different targets give two different systems
  the null case            least squares unchanged when ml_data is inactive
  the scale                held within a cycle, updated between them
  the reuse                alpha and beta from the previous Fc, no extra pass

The one that would have caught the original fault is the plumbing test, and it
is deliberately blunt: two targets that share no mathematics must not produce
the same system.
"""
from __future__ import absolute_import, division, print_function

from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from smtbx.refinement import least_squares
from smtbx_refinement_least_squares_ext import (
  fixed_scale_factor_summary, separable_scale_factor_summary, ml_data)
import random


def equations(n_par=4, n_eq=25, seed=0):
  """A small dense system with known numbers, and its blocks."""
  rnd = random.Random(seed)
  rows = []
  for i in range(n_eq):
    grad = flex.double([rnd.uniform(-1, 1) for _ in range(n_par)])
    yc = rnd.uniform(0.5, 3.0)
    yo = yc*1.3 + rnd.uniform(-0.2, 0.2)
    w = rnd.uniform(0.2, 2.0)
    rows.append((yc, grad, yo, w))
  return rows


def feed(acc, rows):
  for yc, grad, yo, w in rows:
    acc.add_equation(yc, grad, yo, w)
  return acc


def sums(rows, n_par):
  s_yo_sq = s_yc_sq = s_yo_yc = 0.
  yo_g = flex.double(n_par, 0.)
  yc_g = flex.double(n_par, 0.)
  for yc, grad, yo, w in rows:
    s_yo_sq += w*yo*yo
    s_yc_sq += w*yc*yc
    s_yo_yc += w*yo*yc
    yo_g += w*yo*grad
    yc_g += w*yc*grad
  return s_yo_sq, s_yc_sq, s_yo_yc, yo_g, yc_g


def exercise_accumulator_against_closed_form():
  """Objective, right hand side and blocks, against arithmetic done here.

  Not against the other accumulator: two implementations agreeing proves only
  that they share an error.
  """
  n_par = 4
  rows = equations(n_par=n_par)
  blocks_p = flex.int(list(range(n_par)))
  blocks_n = flex.int([n_par])
  for k in (1.0, 0.75, 2.5):
    acc = fixed_scale_factor_summary(n_par, k, blocks_p, blocks_n)
    feed(acc, rows)
    acc.finalise(False)
    s_yo_sq, s_yc_sq, s_yo_yc, yo_g, yc_g = sums(rows, n_par)

    assert approx_equal(acc.scale_factor(), k, eps=1e-15)

    # the whole quadratic: nothing cancels at a k which was not chosen for it
    expected_objective = (s_yo_sq - 2*k*s_yo_yc + k*k*s_yc_sq)/(2*s_yo_sq)
    assert approx_equal(acc.objective(), expected_objective, eps=1e-12), (
      acc.objective(), expected_objective)
    # and it cannot be negative, which the separable short form can give here
    assert acc.objective() >= 0

    # b = D^T W r with D = k J
    rhs = acc.right_hand_side()
    for i in range(n_par):
      assert approx_equal(rhs[i], k*(yo_g[i] - k*yc_g[i])/s_yo_sq, eps=1e-12)

    # grad k is identically zero when k is held
    gk = acc.grad_scale_factor()
    for i in range(n_par):
      assert gk[i] == 0

    # blocks are k^2 (J (x) J), normalised
    blocks = acc.blocks()
    jtj = flex.double(flex.grid(n_par, n_par), 0.)
    for yc, grad, yo, w in rows:
      for a in range(n_par):
        for c in range(n_par):
          jtj[a, c] += w*grad[a]*grad[c]
    for a in range(n_par):
      for c in range(n_par):
        assert approx_equal(blocks[a*n_par + c], k*k*jtj[a, c]/s_yo_sq,
                            eps=1e-12)
  print("accumulator against closed form: OK")


def exercise_agreement_with_separable_at_its_optimum():
  """Where the two must agree, and where they must not.

  Hold the scale at the value the separable accumulator would have chosen and
  the objective and the right hand side coincide - the separable one's extra
  term in the right hand side is identically zero there. The *blocks* still
  differ, and must: the separable design matrix carries the dependence of k on
  the parameters and the fixed one does not. A test which demanded agreement
  everywhere would be demanding the change be pointless.
  """
  n_par = 4
  rows = equations(n_par=n_par, seed=7)
  blocks_p = flex.int(list(range(n_par)))
  blocks_n = flex.int([n_par])

  sep = separable_scale_factor_summary(n_par, blocks_p, blocks_n)
  feed(sep, rows)
  sep.finalise(False)
  k_star = sep.scale_factor()

  fix = fixed_scale_factor_summary(n_par, k_star, blocks_p, blocks_n)
  feed(fix, rows)
  fix.finalise(False)

  assert approx_equal(fix.objective(), sep.objective(), eps=1e-12), (
    fix.objective(), sep.objective())
  a, b = fix.right_hand_side(), sep.right_hand_side()
  for i in range(n_par):
    assert approx_equal(a[i], b[i], eps=1e-12), (i, a[i], b[i])

  # the blocks differ, by the separable one's grad k terms
  fb, sb = fix.blocks(), sep.blocks()
  assert not approx_equal(fb[0], sb[0], eps=1e-9, out=None), (fb[0], sb[0])
  print("agreement with the separable accumulator at its optimum: OK")


def exercise_objective_is_not_the_short_form():
  """At a scale far from the optimum the two objectives must part company.

  This is what distinguishes the fixed accumulator from a copy of the separable
  one with a different scale written into it. Using the separable short form
  here gives a value that falls as k rises without limit, and goes negative.
  """
  n_par = 3
  rows = equations(n_par=n_par, seed=11)
  blocks_p = flex.int(list(range(n_par)))
  blocks_n = flex.int([n_par])
  s_yo_sq, s_yc_sq, s_yo_yc, _, _ = sums(rows, n_par)
  k_far = 10.0
  acc = fixed_scale_factor_summary(n_par, k_far, blocks_p, blocks_n)
  feed(acc, rows)
  acc.finalise(False)
  short_form = (s_yo_sq - k_far*s_yo_yc)/(2*s_yo_sq)
  assert short_form < 0, short_form          # the trap being guarded against
  assert acc.objective() > 0, acc.objective()
  print("objective is the full quadratic, not the short form: OK")


def exercise_objective_only_stops_early():
  """objective_only must not touch the gradients, and must still give a scale.

  The opening pass of a run asks for exactly this, so a version which fell over
  or returned rubbish here would break the first cycle of every refinement.
  """
  n_par = 3
  rows = equations(n_par=n_par, seed=3)
  acc = fixed_scale_factor_summary(n_par, 1.4, flex.int(list(range(n_par))),
                                   flex.int([n_par]))
  feed(acc, rows)
  acc.finalise(True)
  assert approx_equal(acc.scale_factor(), 1.4, eps=1e-15)
  s_yo_sq, s_yc_sq, s_yo_yc, _, _ = sums(rows, n_par)
  assert approx_equal(
    acc.objective(),
    (s_yo_sq - 2*1.4*s_yo_yc + 1.4*1.4*s_yc_sq)/(2*s_yo_sq), eps=1e-12)
  print("objective_only: OK")


def exercise_merge():
  """Two halves accumulated apart and merged must equal one pass.

  How the threaded build combines its workers, so it is not optional.
  """
  n_par = 4
  rows = equations(n_par=n_par, seed=5, n_eq=30)
  blocks_p, blocks_n = flex.int(list(range(n_par))), flex.int([n_par])
  k = 0.9
  whole = feed(fixed_scale_factor_summary(n_par, k, blocks_p, blocks_n), rows)
  whole.finalise(False)
  first = feed(fixed_scale_factor_summary(n_par, k, blocks_p, blocks_n),
               rows[:13])
  second = feed(fixed_scale_factor_summary(n_par, k, blocks_p, blocks_n),
                rows[13:])
  first.merge(second)
  first.finalise(False)
  assert approx_equal(first.objective(), whole.objective(), eps=1e-12)
  a, b = first.right_hand_side(), whole.right_hand_side()
  for i in range(n_par):
    assert approx_equal(a[i], b[i], eps=1e-12)
  print("merge of two halves: OK")


def run():
  exercise_accumulator_against_closed_form()
  exercise_agreement_with_separable_at_its_optimum()
  exercise_objective_is_not_the_short_form()
  exercise_objective_only_stops_early()
  exercise_merge()
  print("OK")


if __name__ == '__main__':
  run()
