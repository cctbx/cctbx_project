from __future__ import absolute_import, division, print_function
import numpy as np
import mmtbx.refinement.llgi_e_dmodel as dmodel
from libtbx.test_utils import approx_equal

def _random_theta_and_grid(rnd, k):
  a = rnd.uniform(0.05, 0.9, size=k)
  b_k_grid = np.sort(rnd.uniform(1.0, 120.0, size=k))
  b = rnd.uniform(0.02, 0.5)
  b_defect = rnd.uniform(5.0, 200.0)
  theta = np.empty(k + 2, dtype=float)
  theta[:k] = a
  theta[-2] = b
  theta[-1] = b_defect
  return theta, b_k_grid

def exercise_unpack_theta_roundtrip():
  rnd = np.random.RandomState(0)
  for k in [0, 1, 2, 3]:
    if(k > 0):
      theta, _ = _random_theta_and_grid(rnd, k)
    else:
      theta = np.array([0.1, 20.0])
    a, b, b_defect = dmodel.unpack_theta(theta)
    assert a.size == k
    if(k > 0):
      assert approx_equal(list(a), list(theta[:k]))
    assert approx_equal(b, theta[-2])
    assert approx_equal(b_defect, theta[-1])

def exercise_unpack_theta_rejects_too_short():
  try:
    dmodel.unpack_theta(np.array([0.1]))
  except ValueError:
    pass
  else:
    raise RuntimeError("expected ValueError for length-1 theta")

def exercise_d_model_raw_zero_a_and_b_gives_zero():
  # With every a_k=0 and b=0, D_raw (the unwrapped, unbounded sum) must
  # be identically zero regardless of the (otherwise-meaningless) fixed
  # b_k_grid/B_defect values -- and D_model = tanh(smooth_relu(D_raw))
  # then follows: smooth_relu(0) is NOT exactly 0 (it is 0.5*sqrt(eps),
  # a deliberately tiny residual -- see the module's own docstring), so
  # this checks D_model is close to (not exactly) 0, well within the
  # documented _SMOOTH_RELU_EPS-driven tolerance, and nowhere near the
  # WRONG 0.5 an earlier, buggy 0.5*(1+tanh(.)) wrapping produced here
  # (caught only after publishing a curve-shape comparison whose every
  # curve plateaued at exactly 0.5 at high resolution -- see design doc
  # sec. 6.4).
  theta = np.array([0.0, 0.0, 0.0, 50.0])  # a_1=a_2=0, b=0, B_defect=50
  b_k_grid = np.array([30.0, 80.0])
  s2 = np.linspace(0.0, 1.0, 25)
  draw = dmodel.d_model_raw(s2, theta, b_k_grid)
  assert flex_max_abs(draw) < 1.e-12
  d = dmodel.d_model(s2, theta, b_k_grid)
  assert flex_max_abs(d) < 1.e-2, d
  assert flex_max_abs(d - 0.5) > 0.4, (
    "D_model at D_raw=0 must be close to 0, not 0.5 -- regression "
    "check for the exact bug found via the published curve comparison")

def exercise_d_model_raw_single_term_matches_hand_formula():
  # b=0 isolates the single coordinate-error term: theta = [a_1, b=0,
  # B_defect] (b/B_defect are mandatory -- see unpack_theta), b_k_grid
  # = [B_1], so D_raw(s) = a_1*exp(-B_1*s2).
  theta = np.array([0.6, 0.0, 1.0])
  b_k_grid = np.array([40.0])
  s2 = np.array([0.0, 0.05, 0.2, 0.5, 1.3])
  draw = dmodel.d_model_raw(s2, theta, b_k_grid)
  expected = 0.6 * np.exp(-40.0 * s2)
  assert flex_max_abs(draw - expected) < 1.e-12
  # D_model = tanh(smooth_relu(D_raw)) -- check the wrapping is
  # applied, not bypassed, on this same case.
  d = dmodel.d_model(s2, theta, b_k_grid)
  eps = 1.e-6  # must match llgi_e_dmodel._SMOOTH_RELU_EPS
  smooth_relu_expected = 0.5*(expected + np.sqrt(expected*expected + eps))
  assert flex_max_abs(d - dmodel.SIGMAA_MAX * np.tanh(smooth_relu_expected)) < 1.e-8

def exercise_d_model_decays_to_zero_at_high_resolution():
  # The defining physical property this module's tanh(smooth_relu(.))
  # wrapping exists to guarantee (design doc sec. 6.4): D_model -> 0,
  # NOT some other constant, as s -> infinity, for ANY theta -- every
  # term in D_raw is a decaying exponential, so D_raw(s->infinity)=0
  # identically, and the wrapping function must send that to 0. A
  # first-attempt fix (0.5*(1+tanh(.))) satisfied non-negativity but
  # sent D_raw=0 to D_model=0.5, so EVERY fitted curve incorrectly
  # plateaued at 0.5 at high resolution -- caught only after the user
  # noticed a published comparison of real fitted curves all showing
  # this exact wrong shape. This test locks the correct asymptote in
  # directly, across a spread of theta (not just the all-zero case
  # above), so this specific regression cannot reappear silently.
  rnd = np.random.RandomState(11)
  s2_huge = np.array([50.0, 200.0, 1.e4])
  for _ in range(5):
    theta, b_k_grid = _random_theta_and_grid(rnd, 2)
    d = dmodel.d_model(s2_huge, theta, b_k_grid)
    assert flex_max_abs(d) < 1.e-2, (theta, b_k_grid, d)

def exercise_d_model_raw_defect_term_matches_hand_formula():
  # a_1=0 isolates the defect term: D_raw(s) = -b*exp(-B_defect*s2).
  theta = np.array([0.0, 0.3, 60.0])
  b_k_grid = np.array([1.0])
  s2 = np.array([0.0, 0.1, 0.4, 0.9])
  draw = dmodel.d_model_raw(s2, theta, b_k_grid)
  expected = -0.3 * np.exp(-60.0 * s2)
  assert flex_max_abs(draw - expected) < 1.e-12

def exercise_d_model_is_bounded_to_open_unit_interval():
  # The defining property of the tanh(smooth_relu(.)) wrapping (design
  # doc sec. 6.4 / this module's own docstring): D_model must be
  # strictly within (0, 1) -- NOT (-1, 1): D_obs and D_model are never
  # negative for any case of practical interest -- for ANY theta,
  # including values that push D_raw far outside any physically sane
  # range (observed in practice during an unconstrained L-BFGS fit's
  # line-search probing on real 2G38 data).
  rnd = np.random.RandomState(9)
  s2 = np.linspace(0.0, 2.0, 15)
  for scale in [1.0, 10.0, 1.e3, 1.e6, 1.e12]:
    theta, b_k_grid = _random_theta_and_grid(rnd, 2)
    theta = theta * scale
    d = dmodel.d_model(s2, theta, b_k_grid)
    assert np.all(np.isfinite(d)), (scale, d)
    assert np.all(d > 0.0), (scale, d)
    assert np.all(d <= dmodel.SIGMAA_MAX), (scale, d)

def exercise_d_model_stays_nonnegative_when_defect_term_dominates():
  # The other defining property of tanh(smooth_relu(.)) (design doc
  # sec. 6.4): D_model must stay non-negative even where the defect
  # term b*exp(-B_defect*s2) exceeds the positive coordinate-error
  # terms at some s (D_raw itself dips below 0 there) -- nothing in
  # the a_k>=0/b>=0/B_defect>=0 constraints alone prevents this, so the
  # wrapping function has to handle it. Plain (unrescaled) tanh(D_raw)
  # would go genuinely negative here; smooth_relu clamps it toward
  # (not exactly to) 0 instead.
  theta = np.array([0.01, 0.01, 5.0, 30.0])  # tiny a_k, huge b
  b_k_grid = np.array([20.0, 50.0])
  s2 = np.array([0.001, 0.01, 0.05])
  draw = dmodel.d_model_raw(s2, theta, b_k_grid)
  assert np.any(draw < 0.0), (
    "test setup should produce a negative D_raw somewhere", draw)
  d = dmodel.d_model(s2, theta, b_k_grid)
  assert np.all(d >= 0.0), d

def exercise_d_model_gradient_matches_finite_difference():
  rnd = np.random.RandomState(1)
  theta, b_k_grid = _random_theta_and_grid(rnd, 2)
  s2 = np.array([0.001, 0.02, 0.1, 0.3, 0.7, 1.5])
  grad = dmodel.d_model_gradient(s2, theta, b_k_grid)
  h = 1.e-6
  for i in range(theta.size):
    tp = theta.copy(); tp[i] += h
    tm = theta.copy(); tm[i] -= h
    fd = (dmodel.d_model(s2, tp, b_k_grid)
          - dmodel.d_model(s2, tm, b_k_grid)) / (2*h)
    rel = np.abs(grad[i] - fd) / np.maximum(1.0, np.abs(fd))
    assert rel.max() < 1.e-5, (i, rel.max())

def exercise_d_model_gradient_is_finite_and_warning_free_for_extreme_theta():
  # Real bug found on real 2G38 data (doc/llgi_target_design.md sec.
  # 6.4): clipping D_raw before tanh() correctly keeps D_model itself
  # bounded (exercise_d_model_is_bounded_to_open_unit_interval), but
  # d_model_raw_gradient is computed from the
  # UNCLIPPED theta and can themselves overflow to inf for extreme
  # theta -- multiplying that against the (correctly tiny, but no
  # longer coupled to the unclipped D_raw) sech2 factor does NOT
  # reproduce the true mathematical limit (which is exactly 0), it
  # produces inf/nan. The fix masks the gradient to exactly 0
  # wherever D_raw was clipped, and suppresses (not hides -- the
  # overflow is provably benign once masked) the resulting numpy
  # RuntimeWarning. This test enforces BOTH: no warning escapes, and
  # the result is finite.
  import warnings
  theta = np.array([1.e300, 0.5, 30.0])  # a_1=1e300, b=0.5, B_defect=30
  b_k_grid = np.array([5.0])
  s2 = np.array([0.001, 0.05, 0.3, 0.7, 1.5])
  with warnings.catch_warnings():
    warnings.simplefilter("error", RuntimeWarning)
    grad = dmodel.d_model_gradient(s2, theta, b_k_grid)
  assert np.all(np.isfinite(grad)), grad
  # D_raw is saturated at every one of these s2 values for a_1=1e300
  # (B_1=5.0 decays far too slowly to bring it back under _D_RAW_CLIP
  # at any of these resolutions), so the a_1 gradient row must be
  # exactly (not approximately) zero.
  assert flex_max_abs(grad[0]) == 0.0, grad[0]

def flex_max_abs(arr):
  return float(np.max(np.abs(arr)))

def run():
  exercise_unpack_theta_roundtrip()
  exercise_unpack_theta_rejects_too_short()
  exercise_d_model_raw_zero_a_and_b_gives_zero()
  exercise_d_model_raw_single_term_matches_hand_formula()
  exercise_d_model_raw_defect_term_matches_hand_formula()
  exercise_d_model_decays_to_zero_at_high_resolution()
  exercise_d_model_stays_nonnegative_when_defect_term_dominates()
  exercise_d_model_is_bounded_to_open_unit_interval()
  exercise_d_model_gradient_matches_finite_difference()
  exercise_d_model_gradient_is_finite_and_warning_free_for_extreme_theta()
  print("OK")

if (__name__ == "__main__"):
  run()
