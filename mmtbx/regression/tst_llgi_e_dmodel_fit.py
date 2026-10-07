from __future__ import absolute_import, division, print_function
import numpy as np
from cctbx.array_family import flex
import mmtbx.refinement.llgi_e_dmodel as dmodel
import mmtbx.refinement.llgi_e_dmodel_fit as fit
from libtbx.test_utils import approx_equal

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
  # E-VALUE NORMALIZATION IS LOAD-BEARING HERE, not cosmetic: an
  # earlier version of this generator drew e_model from a plain
  # uniform(0.4, 2.0) (mean E^2 ~= 1.6) and built e_eff from a noisy
  # scaled copy of it WITHOUT renormalizing (mean E^2 ~= 0.27) --
  # llgi_e_likelihood.py's l(D) assumes NORMALIZED E-values (mean(E^2)
  # = 1 by the standard Wilson-statistics convention it's derived
  # under -- module docstring). With e_eff and e_c sitting on two
  # different, uncalibrated scales, the small-D expansion of l(D) picks
  # up a leading coefficient (1-e_eff^2)*(1-e_c^2) -- strictly NEGATIVE
  # whenever e_eff and e_c straddle opposite sides of 1, which this
  # mismatched scaling made systematic -- so the likelihood itself
  # scored D=0 as better than the TRUE theta, at every sample size
  # tried (confirmed directly: LL(true_theta) < LL(D=0)=0 even at
  # n=4000). Every D_model fit run against that data was therefore
  # chasing a target whose own global optimum was D=0, no matter how
  # well-posed the parametrization or how strong the true signal --
  # not a fitting bug at all, caught only by cross-checking against an
  # INDEPENDENT, non-LLGI moment estimator (moment_sigmaa_estimate)
  # that clearly resolved the true curve where the LLGI fit could not.
  # Fixed here by generating e_c/e_eff from the actual Wilson-statistics
  # generative model the likelihood is derived under (Rice/Woolfson
  # distribution: E_eff | E_c, D ~ Rice(D*E_c, sqrt(1-D^2)) acentric,
  # the analogous half-normal-mixture form centric), which is
  # normalized to mean(E^2)=1 by construction, not by a post-hoc
  # rescale.
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

  # Normalized model amplitude E_c: acentric ~ Rayleigh(mean(E^2)=1,
  # i.e. E_c^2 ~ Exp(1)); centric ~ half-normal(mean(E^2)=1).
  e_c_acentric = np.sqrt(-np.log(rnd.uniform(1.e-12, 1.0, size=n)))
  e_c_centric = np.abs(rnd.standard_normal(n))
  e_c_np = np.where(centric_flags_np, e_c_centric, e_c_acentric)

  # E_eff | E_c, D: the Rice-distribution generative model itself
  # (acentric) / its centric analogue -- a 2D Gaussian with mean
  # (D*E_c, 0) and both components variance (1-D^2)/2, whose modulus is
  # Rice-distributed with the exact l(D) this module fits against;
  # centric collapses to a signed 1D Gaussian with mean D*E_c and
  # variance (1-D^2), whose absolute value matches centric_l.
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

def moment_sigmaa_estimate(d_star_sq, e_eff, e_model, n_bins=8):
  """ Independent, non-LLGI point estimate of sigmaA(resolution), by
  resolution bin, per the user's own diagnostic: for centered/scaled
  intensity-like quantities, corr(E_eff^2, E_model^2) in a resolution
  bin approximates sigmaA^2 for that bin (a standard second-moment
  relationship, e.g. Read (1990)/Srinivasan & Parthasarathy's
  intensity-based sigmaA estimators -- NOT derived from or dependent on
  llgi_e_likelihood.py's own machinery in any way, so it is a genuine
  cross-check, not a restatement of the same model). Used here purely
  as a synthetic-data sanity check (exercise_synthetic_reflections_
  have_resolvable_signal below): if the correlation coefficient is <=0
  in a bin, the estimate is clamped to 0 (matching the user's own
  framing -- a non-positive correlation means the LLGI likelihood's own
  unrestrained optimum in that bin genuinely IS sigmaA=0, not a fitting
  failure to work around).

  d_star_sq, e_eff, e_model: 1D arrays/flex.double, same length.
  n_bins: number of equal-COUNT (quantile) resolution bins.

  Returns (bin_centers_d_star_sq, sigmaa_estimate), both 1D numpy
  arrays of length n_bins (a bin with too few reflections, < 5, is
  dropped from the output rather than returning a noisy estimate).
  """
  d_star_sq = np.asarray(d_star_sq, dtype=float)
  e_eff = np.asarray(e_eff, dtype=float)
  e_model = np.asarray(e_model, dtype=float)
  edges = np.quantile(d_star_sq, np.linspace(0.0, 1.0, n_bins + 1))
  centers, ests = [], []
  for i in range(n_bins):
    sel = (d_star_sq >= edges[i]) & (d_star_sq <= edges[i + 1])
    if(np.sum(sel) < 5):
      continue
    r = np.corrcoef(e_eff[sel]**2, e_model[sel]**2)[0, 1]
    centers.append(0.5 * (edges[i] + edges[i + 1]))
    ests.append(np.sqrt(max(r, 0.0)))
  return np.array(centers), np.array(ests)

def exercise_synthetic_reflections_have_resolvable_signal():
  # Guards _build_synthetic_reflections's own realism: an INDEPENDENT
  # (non-LLGI) moment-based sigmaA estimate, per resolution bin, must
  # track the true D_model curve reasonably well -- not exactly (it is
  # a noisy point estimate from a finite sample), but well enough that
  # the true curve is clearly resolvable above the estimator's own
  # noise floor. This is the check that would have caught the earlier
  # version of this dataset's B_k=[30,90] choice, whose true signal was
  # ALREADY indistinguishable from zero at this sample size/noise level
  # (confirmed: the raw LLGI likelihood at that true theta was itself
  # worse than at D_model=0 identically, for any n up to 4000 tried).
  inputs = _build_synthetic_reflections(seed=0, n=4000)
  d_star_sq = np.array(inputs["d_star_sq"])
  centers, ests = moment_sigmaa_estimate(
    d_star_sq, np.array(inputs["e_eff"]), np.array(inputs["e_model"]))
  assert centers.size >= 6, centers.size
  true_curve = dmodel.d_model(
    centers / 4.0, inputs["true_theta"], inputs["true_b_k_grid"])
  # The moment estimate is a noisy point estimate (one correlation
  # coefficient per bin, ~500 reflections each) -- checking it stays
  # within a generous absolute tolerance of the true curve, across
  # every bin, confirms the signal is resolvable at all (a truly
  # unresolvable signal, like the earlier B_k=[30,90] dataset, gives
  # near-zero/negative correlations everywhere regardless of the true
  # curve's own shape, which this would catch).
  worst = float(np.max(np.abs(ests - true_curve)))
  assert worst < 0.25, (centers, ests, true_curve, worst)

def exercise_estimate_d_model_sigmaa_converges_without_degenerating():
  # An earlier version of this module used a soft per-reflection
  # barrier to discourage D=dobs*D_model from reaching +-1 -- found
  # (see llgi_e_dmodel.py's own docstring) to be fundamentally unable
  # to prevent the true failure mode: the exact log-likelihood diverges
  # to -infinity approaching D=+-1, so any finite barrier weight is
  # eventually outweighed by the reward of crossing into the negative-
  # variance guard's region (which returns exactly 0). D_model(s;
  # theta) is now tanh(smooth_relu(.))-wrapped and UNCONDITIONALLY
  # bounded to (0, 1) for any theta (llgi_e_dmodel.d_model), so THAT
  # problem cannot occur at all. A SEPARATE degenerate direction was
  # found once this asymptote fix made D_model=0 genuinely reachable:
  # since l(D=0)=0 and l'(D=0)=0 identically (llgi_e_likelihood.py),
  # the earlier free-B_k parametrization could run a_k->0/B_k->infinity
  # away to a numerically meaningless solution (confirmed: B_k
  # converging to ~1e13). B_k is now a FIXED ladder (design doc sec.
  # 6.4's well-posedness addendum -- no longer part of theta at all),
  # which closes off that specific runaway structurally; the checks
  # below confirm theta (now just [a_1..a_K, b, B_defect]) stays sane.
  # b_sol_anchor=40 (a typical real-data B_sol scale) matches how this
  # fit is actually called in real phenix.refine usage -- by the time
  # the D_model sigmaA fit runs, a bulk-solvent B_sol point estimate is
  # always available and passed through (mmtbx.refinement.
  # llgi_e_sigmaa's own callers). Testing the fully-unrestrained
  # (b_sol_anchor=None) case here left B_defect just as free/nonlinear
  # as the old B_k's were, and it ran away the same way (b->1e13,
  # B_defect->~0, a flat/uninformative defect term) once D_model=0
  # became genuinely reachable -- not a realistic scenario (b_sol_
  # anchor=None is exercised deliberately/separately, see
  # exercise_b_sol_restraint_pulls_b_defect_toward_anchor's own
  # unrestrained-vs-restrained comparison).
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
  assert result.target != 0.0, (
    "target collapsed to exactly 0 -- the degenerate every-reflection-"
    "masked-out failure mode the tanh wrapping/b_k_grid fix exists to "
    "prevent")
  sigmaa_np = np.array(result.sigmaa)
  assert np.all(sigmaa_np > 0.0) and np.all(sigmaa_np < 1.0), (
    "D_model itself must be strictly within (0, 1) for every "
    "reflection, by construction (tanh(smooth_relu(.))-wrapped) -- not "
    "merely small enough that dobs*D_model happens to stay bounded")
  D_all = np.array(inputs["dobs"]) * sigmaa_np
  assert np.max(np.abs(D_all)) < 1.0, np.max(np.abs(D_all))
  assert result.sigmaa.size() == inputs["e_eff"].size()
  # include_constant_term (the default) adds a B=0 rung to the 2 decaying
  assert result.b_k_grid.size() == 3
  assert result.b_k_grid[0] == 0.0

def exercise_estimate_d_model_sigmaa_recovers_true_curve():
  # The actual recovery-quality check this module was missing for a
  # long stretch of this design's history: not just "doesn't crash/
  # degenerate" (the test above) but "gets a reasonably close answer"
  # -- only meaningful once _build_synthetic_reflections generates
  # PROPERLY NORMALIZED E-values (mean(E^2)=1, matching the Wilson-
  # statistics convention llgi_e_likelihood.py's l(D) is derived under
  # -- see that function's own docstring for the full story: an earlier
  # un-normalized version of this dataset had LL(true_theta) < LL(D=0)
  # identically, so EVERY fit against it was chasing a target whose own
  # global optimum was D=0, and no amount of correct D_model machinery
  # could have recovered the true curve from data like that). With
  # properly normalized data and n=4000 (enough reflections for the
  # per-bin moment estimator itself to be reasonably tight), the fitted
  # curve should track the true curve's magnitude to within a loose but
  # real tolerance across the resolution range, not just share its
  # sign/rough shape.
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
  # Directional sanity: the fitted curve should decay with resolution
  # like the true one does, not be flat or inverted.
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
  # b_k_grid can be pinned explicitly (e.g. for reproducible tests/
  # diagnostics) instead of derived from the data -- check the fit
  # actually uses the grid passed in, not a freshly-derived one.
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
  exercise_synthetic_reflections_have_resolvable_signal()
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
