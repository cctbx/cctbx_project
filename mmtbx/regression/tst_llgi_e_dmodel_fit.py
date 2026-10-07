from __future__ import absolute_import, division, print_function
import numpy as np
from cctbx.array_family import flex
import mmtbx.refinement.llgi_e_dmodel as dmodel
import mmtbx.refinement.llgi_e_dmodel_fit as fit
from libtbx.test_utils import approx_equal
from mmtbx.regression.llgi_test_utils import random_theta_and_grid

def exercise_default_b_k_grid_spans_data_resolution_range():
  # default_b_k_grid picks a log-spaced ladder from roughly 1/s2_max to
  # 1/s2_min (llgi_e_dmodel_fit.default_b_k_grid's own docstring: the
  # B_k range at which exp(-B_k*s2) actually varies somewhere within
  # the data). Check the ladder brackets that range and is monotonic.
  rnd = np.random.RandomState(5)
  s2 = rnd.uniform(0.01, 0.8, size=500)
  grid = fit.default_b_k_grid(4, s2)
  assert grid.shape == (4,)
  assert np.all(np.diff(grid) > 0.0), grid  # strictly ascending
  b_lo_expected = 1.0 / s2.max()
  b_hi_expected = 1.0 / s2.min()
  assert grid[0] >= b_lo_expected * 0.5, (grid, b_lo_expected)
  assert grid[-1] <= b_hi_expected * 2.0, (grid, b_hi_expected)

def exercise_default_b_k_grid_falls_back_for_degenerate_s2():
  grid_empty = fit.default_b_k_grid(3, np.array([]))
  assert grid_empty.shape == (3,)
  assert np.all(np.isfinite(grid_empty))
  grid_constant = fit.default_b_k_grid(3, np.array([0.5, 0.5, 0.5]))
  assert grid_constant.shape == (3,)
  assert np.all(np.isfinite(grid_constant))

def exercise_default_b_k_grid_k_zero_and_one():
  assert fit.default_b_k_grid(0, np.array([0.1, 0.5])).shape == (0,)
  grid_one = fit.default_b_k_grid(1, np.array([0.1, 0.5]))
  assert grid_one.shape == (1,)
  assert grid_one[0] > 0.0

def _build_synthetic_reflections(seed, n=400, true_theta=None,
      true_b_k_grid=None):
  # Normalised E values from the model the likelihood assumes:
  # |E_c|^2 ~ Exp(1) (acentric) or half-normal E_c (centric), and E_eff
  # Rice-distributed about D*E_c with variance 1 - D^2 (acentric), or
  # |N(D*E_c, 1 - D^2)| (centric).
  rnd = np.random.RandomState(seed)
  d_star_sq = flex.double(rnd.uniform(0.02, 3.0, size=n))
  s2 = np.array(d_star_sq) / 4.0
  if(true_b_k_grid is None):
    true_b_k_grid = np.array([3.0, 15.0])
  if(true_theta is None):
    true_theta = np.array([0.7, 0.5, 0.1, 20.0])  # a_1, a_2, b, B_defect
  true_b_k_grid = np.asarray(true_b_k_grid, dtype=float)
  true_theta = np.asarray(true_theta, dtype=float)
  true_dmodel = dmodel.d_model(s2, true_theta, true_b_k_grid)
  dobs_true_np = np.clip(0.85 + 0.1*rnd.standard_normal(n), 0.3, 0.98)
  D_true = dobs_true_np * true_dmodel
  centric_flags_np = rnd.uniform(size=n) < 0.2

  e_c_acentric = np.sqrt(-np.log(rnd.uniform(1.e-12, 1.0, size=n)))
  e_c_centric = np.abs(rnd.standard_normal(n))
  e_c_np = np.where(centric_flags_np, e_c_centric, e_c_acentric)

  var = np.clip(1.0 - D_true*D_true, 1.e-6, None)
  real = D_true*e_c_np + np.sqrt(var/2.0)*rnd.standard_normal(n)
  imag = np.sqrt(var/2.0)*rnd.standard_normal(n)
  e_eff_acentric = np.sqrt(real*real + imag*imag)
  e_eff_centric = np.abs(D_true*e_c_np + np.sqrt(var)*rnd.standard_normal(n))
  e_eff_np = np.where(centric_flags_np, e_eff_centric, e_eff_acentric)

  e_c = flex.double(e_c_np)
  dobs_true = flex.double(dobs_true_np)
  e_eff = flex.double(e_eff_np)
  r_free_flags = flex.bool(rnd.uniform(size=n) < 0.5)
  centric_flags = flex.bool(centric_flags_np)
  return dict(
    e_eff=e_eff, r_free_flags=r_free_flags, e_model=e_c, dobs=dobs_true,
    centric_flags=centric_flags, d_star_sq=d_star_sq,
    true_theta=true_theta, true_b_k_grid=true_b_k_grid)

def exercise_target_and_gradient_finite_difference():
  # Mixed acentric/centric reflections, no hybrid.
  rnd = np.random.RandomState(11)
  n = 30
  s2 = rnd.uniform(0.0, 1.2, size=n)
  args = (s2, flex.double(rnd.uniform(0.3, 2.2, size=n)),
    flex.double(rnd.uniform(0.3, 2.2, size=n)),
    flex.double(rnd.uniform(0.3, 0.95, size=n)),
    flex.bool((rnd.uniform(size=n) < 0.25).tolist()))
  theta, b_k_grid = random_theta_and_grid(rnd, 2)
  t, grad = fit.target_and_gradient(theta, *(args + (b_k_grid,)))
  h = 1.e-6
  for i in range(theta.size):
    tp = theta.copy(); tp[i] += h
    tm = theta.copy(); tm[i] -= h
    fd = (fit.target_and_gradient(tp, *(args + (b_k_grid,)))[0]
        - fit.target_and_gradient(tm, *(args + (b_k_grid,)))[0]) / (2*h)
    assert abs(grad[i] - fd) < 1.e-6 * max(1.0, abs(fd)), (i, grad[i], fd)

def exercise_estimate_d_model_sigmaa_converges_without_degenerating():
  # theta stays finite and within bounds, D_model in (0, 1), and the fit
  # reports normal termination. b_sol_anchor = 40 as in refinement, where
  # bss's B_sol is always passed.
  inputs = _build_synthetic_reflections(seed=0)
  result = fit.estimate_d_model_sigmaa(
    e_eff=inputs["e_eff"], r_free_flags=inputs["r_free_flags"],
    e_model=inputs["e_model"], dobs=inputs["dobs"],
    centric_flags=inputs["centric_flags"],
    d_star_sq=inputs["d_star_sq"], n_gaussian_terms=2, max_iterations=200,
    b_sol_anchor=40.0)
  theta = np.array(result.theta)
  assert np.all(np.isfinite(theta)), theta
  # amplitudes within the fitter's bounds (0 is allowed: a term can drop out)
  assert np.all(theta[:-1] >= 0.0), theta
  assert np.all(theta[:-1] <= fit.d_model_target_evaluator.amplitude_max), theta
  assert result.target < 0.0, result.target
  sigmaa_np = np.array(result.sigmaa)
  assert np.all(sigmaa_np > 0.0) and np.all(sigmaa_np < 1.0)
  D_all = np.array(inputs["dobs"]) * sigmaa_np
  assert np.max(np.abs(D_all)) < 1.0, np.max(np.abs(D_all))
  assert result.lbfgs_error is None, result.lbfgs_error
  assert result.sigmaa.size() == inputs["e_eff"].size()
  # include_constant_term (the default) adds a B=0 rung to the 2 decaying
  assert result.b_k_grid.size() == 3
  assert result.b_k_grid[0] == 0.0

def exercise_estimate_d_model_sigmaa_recovers_true_curve():
  # With n = 4000 the fitted curve follows the true one within 0.3 and
  # decays with resolution.
  inputs = _build_synthetic_reflections(seed=1, n=4000)
  result = fit.estimate_d_model_sigmaa(
    e_eff=inputs["e_eff"], r_free_flags=inputs["r_free_flags"],
    e_model=inputs["e_model"], dobs=inputs["dobs"],
    centric_flags=inputs["centric_flags"],
    d_star_sq=inputs["d_star_sq"], n_gaussian_terms=2, max_iterations=200,
    b_sol_anchor=40.0)
  theta = np.array(result.theta)
  b_k_grid = np.array(result.b_k_grid)
  s2_check = np.array([0.02, 0.1, 0.3, 0.5, 0.75])
  fitted = dmodel.d_model(s2_check, theta, b_k_grid)
  true = dmodel.d_model(s2_check, inputs["true_theta"], inputs["true_b_k_grid"])
  worst = float(np.max(np.abs(fitted - true)))
  assert worst < 0.3, (s2_check, fitted, true, worst)
  assert fitted[0] > fitted[-1], fitted

def exercise_estimate_d_model_sigmaa_no_test_set_raises():
  inputs = _build_synthetic_reflections(seed=2, n=20)
  all_false = flex.bool(inputs["r_free_flags"].size(), False)
  try:
    fit.estimate_d_model_sigmaa(
      e_eff=inputs["e_eff"], r_free_flags=all_false,
      e_model=inputs["e_model"], dobs=inputs["dobs"],
      centric_flags=inputs["centric_flags"],
      d_star_sq=inputs["d_star_sq"], n_gaussian_terms=2)
  except RuntimeError:
    pass
  else:
    raise RuntimeError("expected RuntimeError for an empty test set")

def exercise_estimate_d_model_sigmaa_accepts_explicit_b_k_grid():
  # An explicit b_k_grid is used as given.
  inputs = _build_synthetic_reflections(seed=4)
  explicit_grid = np.array([25.0, 75.0])
  result = fit.estimate_d_model_sigmaa(
    e_eff=inputs["e_eff"], r_free_flags=inputs["r_free_flags"],
    e_model=inputs["e_model"], dobs=inputs["dobs"],
    centric_flags=inputs["centric_flags"],
    d_star_sq=inputs["d_star_sq"], n_gaussian_terms=2, max_iterations=50,
    b_k_grid=explicit_grid, include_constant_term=False)
  assert approx_equal(list(result.b_k_grid), list(explicit_grid))

def exercise_constant_term_ladder_and_fit():
  # include_constant_term prepends a B=0 rung to the data-derived ladder,
  # and on data whose true sigmaA has a resolution-independent component
  # (flat at high resolution, as for a well-refined model) it fits better
  # than the same ladder without it.
  inputs = _build_synthetic_reflections(seed=6, n=4000,
    true_theta=np.array([1.2, 0.6, 0.02, 20.0]),
    true_b_k_grid=np.array([0.0, 10.0]))
  common = dict(
    e_eff=inputs["e_eff"], r_free_flags=inputs["r_free_flags"],
    e_model=inputs["e_model"], dobs=inputs["dobs"],
    centric_flags=inputs["centric_flags"], d_star_sq=inputs["d_star_sq"],
    n_gaussian_terms=2, max_iterations=200, b_sol_anchor=20.0)
  plain = fit.estimate_d_model_sigmaa(include_constant_term=False, **common)
  with_const = fit.estimate_d_model_sigmaa(include_constant_term=True,
    **common)
  grid = np.array(with_const.b_k_grid)
  assert grid.size == 3 and grid[0] == 0.0, grid
  assert approx_equal(list(grid[1:]), list(plain.b_k_grid))
  assert with_const.theta.size() == 3 + 2
  assert with_const.target < plain.target - 1.e-3, (
    with_const.target, plain.target)

def exercise_constant_term_nests_plain_model():
  # The constant-term model contains plain D_model (a_0 = 0), so it must
  # fit at least as well even when the true curve has no constant
  # component.
  inputs = _build_synthetic_reflections(seed=6, n=4000,
    true_theta=np.array([1.2, 0.0, 0.02, 20.0]),
    true_b_k_grid=np.array([10.0, 60.0]))
  common = dict(
    e_eff=inputs["e_eff"], r_free_flags=inputs["r_free_flags"],
    e_model=inputs["e_model"], dobs=inputs["dobs"],
    centric_flags=inputs["centric_flags"], d_star_sq=inputs["d_star_sq"],
    n_gaussian_terms=2, max_iterations=200, b_sol_anchor=20.0)
  plain = fit.estimate_d_model_sigmaa(include_constant_term=False, **common)
  with_const = fit.estimate_d_model_sigmaa(include_constant_term=True,
    **common)
  assert with_const.target <= plain.target + 1.e-5, (
    with_const.target, plain.target)

def exercise_b_defect_fixed_to_anchor():
  # With an anchor, B_defect is fixed at it exactly (not fitted); without
  # one it is fitted.
  inputs = _build_synthetic_reflections(seed=3)
  common = dict(
    e_eff=inputs["e_eff"], r_free_flags=inputs["r_free_flags"],
    e_model=inputs["e_model"], dobs=inputs["dobs"],
    centric_flags=inputs["centric_flags"], d_star_sq=inputs["d_star_sq"],
    n_gaussian_terms=2, max_iterations=200, include_constant_term=False)
  for anchor in (25.0, 200.0):
    result = fit.estimate_d_model_sigmaa(b_sol_anchor=anchor, **common)
    assert result.theta.size() == 2 + 2
    assert result.theta[-1] == anchor, (result.theta[-1], anchor)
  start = fit.d_model_target_evaluator._default_theta_start(2)
  free = fit.estimate_d_model_sigmaa(b_sol_anchor=None, **common)
  assert free.theta.size() == 2 + 2
  assert abs(free.theta[-1] - start[-1]) > 1.e-3, (free.theta[-1], start[-1])

def run():
  exercise_default_b_k_grid_spans_data_resolution_range()
  exercise_default_b_k_grid_falls_back_for_degenerate_s2()
  exercise_default_b_k_grid_k_zero_and_one()
  exercise_target_and_gradient_finite_difference()
  exercise_estimate_d_model_sigmaa_converges_without_degenerating()
  exercise_estimate_d_model_sigmaa_recovers_true_curve()
  exercise_estimate_d_model_sigmaa_no_test_set_raises()
  exercise_estimate_d_model_sigmaa_accepts_explicit_b_k_grid()
  exercise_constant_term_ladder_and_fit()
  exercise_constant_term_nests_plain_model()
  exercise_b_defect_fixed_to_anchor()
  print("OK")

if (__name__ == "__main__"):
  run()
