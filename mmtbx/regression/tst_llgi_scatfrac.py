from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import mmtbx.refinement.llgi_scatfrac as llgi_scatfrac
from libtbx.test_utils import approx_equal
import math
import random

def rice_sample(ec, v, random_state):
  """ One sample from the acentric Rice distribution with location ec and
  total variance v: the modulus of (ec, 0) plus 2D Gaussian noise of
  variance v/2 per component, as llgi.h's target_one_h assumes. """
  sigma = math.sqrt(max(v, 1.e-10) / 2.0)
  re = random_state.gauss(ec, sigma)
  im = random_state.gauss(0.0, sigma)
  return math.sqrt(re * re + im * im)

def synthetic_dataset(d_star_sq, sigmaa_true, scatfrac_true, dobs,
      random_state, fcalc_mean_sq=None):
  """ Data consistent with the F-scale LLGI model (TEPS = RESN = k = 1):
  |Fcalc|^2 exponential with mean ScatFrac (or fcalc_mean_sq), Feff
  Rice-distributed about d*|Fcalc|, d = Dobs*sigmaA/sqrt(ScatFrac), with
  variance 1 - (Dobs*sigmaA)^2. """
  n_refl = d_star_sq.size()
  if(fcalc_mean_sq is None): fcalc_mean_sq = scatfrac_true
  assert flex.max(dobs * sigmaa_true) < 1
  f_calc = flex.complex_double()
  f_eff = flex.double()
  for i in range(n_refl):
    fc_mag = math.sqrt(random_state.expovariate(1.0 / fcalc_mean_sq[i]))
    phase = random_state.uniform(0.0, 2.0 * math.pi)
    f_calc.append(complex(fc_mag * math.cos(phase), fc_mag * math.sin(phase)))
    big_d = dobs[i] * sigmaa_true[i]
    d = big_d / math.sqrt(scatfrac_true[i])
    f_eff.append(rice_sample(d * fc_mag, 1.0 - big_d * big_d, random_state))
  return dict(
    d_star_sq=d_star_sq, sigmaa_true=sigmaa_true,
    scatfrac_true=scatfrac_true, dobs=dobs,
    teps=flex.double(n_refl, 1.0), resn=flex.double(n_refl, 1.0),
    f_calc=f_calc, f_eff=f_eff,
    centric_flags=flex.bool(n_refl, False),
    working_selection=flex.bool(n_refl, True))

def random_d_star_sq(n_refl, random_state, d_lo=0.001, d_hi=0.25):
  return flex.double(sorted(
    random_state.uniform(d_lo, d_hi) for i in range(n_refl)))

def b_factor_dataset(n_refl, seed, scatfrac_inf, b_scatfrac, sigmaa=0.6,
      dobs=0.85):
  """ ScatFrac = scatfrac_inf*exp(-b_scatfrac*d*^2/4), constant sigmaA. """
  random_state = random.Random(seed)
  d_star_sq = random_d_star_sq(n_refl, random_state)
  scatfrac_true = scatfrac_inf * flex.exp(-b_scatfrac * d_star_sq / 4)
  return synthetic_dataset(d_star_sq, flex.double(n_refl, sigmaa),
    scatfrac_true, flex.double(n_refl, dobs), random_state)

def fit(data, params=None, sigmaa=None, dobs=None, hybrid=None):
  return llgi_scatfrac.estimate_llgi_scatfrac_likelihood(
    f_eff=data["f_eff"], working_selection=data["working_selection"],
    f_calc=data["f_calc"],
    dobs=data["dobs"] if dobs is None else dobs,
    sigmaa=data["sigmaa_true"] if sigmaa is None else sigmaa,
    teps=data["teps"], resn=data["resn"],
    centric_flags=data["centric_flags"], d_star_sq=data["d_star_sq"],
    params=params, hybrid=hybrid)

def exercise_ratio_of_sums_robust_to_single_outlier_reflection():
  # RESN varying ~90x within one resolution range (as in real data) and
  # one reflection at the smallest RESN with |Fcalc| 10x its expected
  # value: the ratio of sums stays close to the truth, a mean of
  # per-reflection ratios does not.
  n_refl = 200
  random_state = random.Random(13)
  teps = flex.double(n_refl, 1.0)
  resn = flex.double([
    math.exp(random_state.uniform(math.log(2.0), math.log(180.0)))
    for i in range(n_refl)])
  true_scatfrac = 0.6
  f_calc = flex.complex_double()
  for i in range(n_refl):
    mean_fc_sq = true_scatfrac * resn[i] ** 2
    fc_mag = math.sqrt(random_state.expovariate(1.0 / mean_fc_sq))
    phase = random_state.uniform(0.0, 2.0 * math.pi)
    f_calc.append(complex(fc_mag * math.cos(phase), fc_mag * math.sin(phase)))
  i_min = flex.min_index(resn)
  f_calc[i_min] = complex(10.0 * math.sqrt(true_scatfrac) * resn[i_min], 0)
  estimate = llgi_scatfrac.scatfrac_ratio_of_sums(
    f_calc=f_calc, teps=teps, resn=resn)
  assert abs(estimate - true_scatfrac) < 0.15, estimate
  mean_of_ratios = flex.mean(flex.norm(f_calc) / (teps * resn * resn))
  assert abs(mean_of_ratios - true_scatfrac) > 0.15, mean_of_ratios
  # scale_factor enters squared
  assert approx_equal(llgi_scatfrac.scatfrac_ratio_of_sums(
    f_calc=f_calc, teps=teps, resn=resn, scale_factor=2.0), 4 * estimate)

def exercise_scatfrac_likelihood_uses_working_set_selection():
  # Two halves with different true ScatFrac (0.4 and 1.1, above 1 on
  # purpose): a fit on each half recovers that half's value.
  n_refl_each = 3000
  random_state = random.Random(56)
  def make_half(scatfrac_level):
    return synthetic_dataset(
      random_d_star_sq(n_refl_each, random_state),
      flex.double(n_refl_each, 0.6), flex.double(n_refl_each, scatfrac_level),
      flex.double(n_refl_each, 0.85), random_state)
  low = make_half(0.4)
  high = make_half(1.1)
  data = {}
  for key in low:
    data[key] = low[key].concatenate(high[key])
  select_low = flex.bool([True] * n_refl_each + [False] * n_refl_each)
  params = llgi_scatfrac.llgi_scatfrac_params.extract()
  data["working_selection"] = select_low
  fit_low = fit(data, params)
  data["working_selection"] = ~select_low
  fit_high = fit(data, params)
  mean_low_on_low = flex.mean(fit_low.scatfrac.select(select_low))
  mean_high_on_high = flex.mean(fit_high.scatfrac.select(~select_low))
  assert abs(mean_low_on_low - 0.4) < 0.15, mean_low_on_low
  assert abs(mean_high_on_high - 1.1) < 0.15, mean_high_on_high
  assert fit_low.lbfgs_error is None and fit_high.lbfgs_error is None

  data["working_selection"] = flex.bool(2 * n_refl_each, False)
  try:
    fit(data)
  except RuntimeError as e:
    assert "working-set" in str(e)
  else:
    raise RuntimeError("Expected RuntimeError for an empty working set.")

def exercise_scatfrac_b_factor_evaluator_gradient_finite_difference():
  # Finite-difference check of the evaluator's gradient, with sigmaA from
  # 0.95 at low resolution to 0.35 at high resolution.
  n_refl = 400
  random_state = random.Random(81)
  d_star_sq = random_d_star_sq(n_refl, random_state)
  sigmaa = 0.95 - 0.6 * d_star_sq / 0.25
  data = synthetic_dataset(d_star_sq, sigmaa,
    0.9 - 0.3 * d_star_sq / 0.25, flex.double(n_refl, 0.85), random_state)
  working_selection = flex.bool([(i % 5 != 0) for i in range(n_refl)])
  evaluator = llgi_scatfrac.llgi_scatfrac_b_factor_target_evaluator(
    f_eff=data["f_eff"], selection=working_selection, f_calc=data["f_calc"],
    dobs=data["dobs"], sigmaa=sigmaa, teps=data["teps"],
    resn=data["resn"], ss=d_star_sq / 4, centric_flags=data["centric_flags"],
    scale_factor=1.0, scatfrac_inf_start=0.7, b_scatfrac_start=15.0,
    restraint_sigma=10.0, max_iterations=0)
  # max_iterations=0 still takes one L-BFGS step: set the point to probe
  evaluator.x = evaluator.x_from(0.7, 15.0)
  x0 = evaluator.x.deep_copy()
  f0, g0 = evaluator.compute_functional_and_gradients()
  eps = 1.e-6
  for i in range(2):
    x = x0.deep_copy(); x[i] += eps
    evaluator.x = x
    fp, _ = evaluator.compute_functional_and_gradients()
    x = x0.deep_copy(); x[i] -= eps
    evaluator.x = x
    fm, _ = evaluator.compute_functional_and_gradients()
    assert approx_equal(g0[i], (fp - fm) / (2 * eps), eps=1.e-6), i

def _assert_b_factor_curve_matches_truth_over_observed_range(
      result, d_star_sq, scatfrac_inf_true, b_scatfrac_true):
  # ScatFrac_inf and B_scatfrac are strongly correlated over a finite
  # d*^2 range, so compare the curve over the observed range, and the
  # sign of B_scatfrac, rather than the two parameters.
  true_vals = scatfrac_inf_true * flex.exp(-b_scatfrac_true * d_star_sq / 4)
  mean_rel_diff = flex.mean(flex.abs(true_vals - result.scatfrac) / true_vals)
  assert mean_rel_diff < 0.3, mean_rel_diff
  log_corr = flex.linear_correlation(
    flex.log(true_vals), flex.log(result.scatfrac)).coefficient()
  assert log_corr > 0.999, log_corr
  assert (result.b_scatfrac > 0) == (b_scatfrac_true > 0), (
    result.b_scatfrac, b_scatfrac_true)

def exercise_scatfrac_b_factor_recovers_falling_trend():
  # Positive B_scatfrac, restraint off.
  data = b_factor_dataset(n_refl=6000, seed=83, scatfrac_inf=0.95,
    b_scatfrac=25.0, sigmaa=0.4)
  params = llgi_scatfrac.llgi_scatfrac_params.extract()
  params.scatfrac_b_factor_restraint_sigma = 0.0
  result = fit(data, params)
  _assert_b_factor_curve_matches_truth_over_observed_range(
    result, data["d_star_sq"], 0.95, 25.0)

def exercise_scatfrac_b_factor_recovers_rising_trend():
  # Negative B_scatfrac (ScatFrac rising toward high resolution, as for a
  # partial model of the best-ordered parts), restraint off.
  data = b_factor_dataset(n_refl=6000, seed=89, scatfrac_inf=0.5,
    b_scatfrac=-25.0)
  params = llgi_scatfrac.llgi_scatfrac_params.extract()
  params.scatfrac_b_factor_restraint_sigma = 0.0
  result = fit(data, params)
  _assert_b_factor_curve_matches_truth_over_observed_range(
    result, data["d_star_sq"], 0.5, -25.0)
  i_low = flex.min_index(data["d_star_sq"])
  i_high = flex.max_index(data["d_star_sq"])
  assert result.scatfrac[i_high] > result.scatfrac[i_low]

def exercise_b_factor_restraint_penalty_finite_difference():
  sigma = 10.0
  for b in (-30.0, -1.0, 0.0, 4.0, 25.0):
    f0, g0 = llgi_scatfrac._b_factor_restraint_penalty_and_gradient(b, sigma)
    eps = 1.e-6
    fp, _ = llgi_scatfrac._b_factor_restraint_penalty_and_gradient(
      b + eps, sigma)
    fm, _ = llgi_scatfrac._b_factor_restraint_penalty_and_gradient(
      b - eps, sigma)
    assert approx_equal(g0, (fp - fm) / (2 * eps), eps=1.e-4)
  assert llgi_scatfrac._b_factor_restraint_penalty_and_gradient(25.0, 0.0) \
    == (0.0, 0.0)

def exercise_scatfrac_b_factor_restraint_pulls_toward_zero():
  # With a small true B_scatfrac, a restraint (sigma = 10) gives a
  # |B_scatfrac| no larger than the unrestrained fit (the default).
  data = b_factor_dataset(n_refl=6000, seed=97, scatfrac_inf=0.8,
    b_scatfrac=3.0)
  params_unrestrained = llgi_scatfrac.llgi_scatfrac_params.extract()
  assert params_unrestrained.scatfrac_b_factor_restraint_sigma == 0.0
  result_unrestrained = fit(data, params_unrestrained)
  params_restrained = llgi_scatfrac.llgi_scatfrac_params.extract()
  params_restrained.scatfrac_b_factor_restraint_sigma = 10.0
  result_restrained = fit(data, params_restrained)
  assert abs(result_restrained.b_scatfrac) < abs(
    result_unrestrained.b_scatfrac) + 1.e-6, (
      result_restrained.b_scatfrac, result_unrestrained.b_scatfrac)

def exercise_likelihood_agrees_with_moments():
  # On data consistent with the model (mean |Fcalc|^2 = ScatFrac), the
  # likelihood fit and the moments estimate of ScatFrac agree, at high and
  # moderate sigmaA.
  n_refl = 3000
  for seed, sigmaa in ((101, 0.95), (103, 0.6)):
    random_state = random.Random(seed)
    d_star_sq = random_d_star_sq(n_refl, random_state)
    scatfrac_true = 0.8 * flex.exp(-6.0 * d_star_sq / 4)
    data = synthetic_dataset(d_star_sq, flex.double(n_refl, sigmaa),
      scatfrac_true, flex.double(n_refl, 0.9), random_state)
    result = fit(data)
    assert result.lbfgs_error is None, result.lbfgs_error
    s_m = llgi_scatfrac.scatfrac_moments(data["f_calc"], data["teps"],
      data["resn"], d_star_sq, n_bins=10)
    rel = flex.mean(flex.abs(result.scatfrac - s_m) / s_m)
    assert rel < 0.1, (seed, rel)
    rel_true = flex.mean(flex.abs(result.scatfrac - scatfrac_true)
      / scatfrac_true)
    assert rel_true < 0.1, (seed, rel_true)

def exercise_scatfrac_fit_with_hybrid_converges():
  # Synthetic fmodel whose sigmaA reaches 0.995: the fit converges, and
  # no point on a coarse grid around it is better.
  from mmtbx.regression.llgi_test_utils import build_llgi_fmodel
  import mmtbx.refinement.llgi_hybrid as llgi_hybrid
  fmodel = build_llgi_fmodel(80, 1.6, seed=4, rice_kappa=0.1)
  llgi_data = fmodel.llgi_data()
  f_obs = fmodel.f_obs()
  k = fmodel.scale_ml_wrapper()
  args = dict(f_eff=llgi_data.feff.data(),
    working_selection=~fmodel.r_free_flags().data(),
    f_calc=fmodel.f_model().data(), dobs=llgi_data.dobs.data(),
    sigmaa=llgi_data.sigmaa.data(), teps=llgi_data.teps.data(),
    resn=llgi_data.resn.data(), centric_flags=f_obs.centric_flags().data(),
    d_star_sq=f_obs.d_star_sq().data(), scale_factor=k)
  result = llgi_scatfrac.estimate_llgi_scatfrac_likelihood(
    hybrid=llgi_hybrid.get_hybrid(llgi_data), **args)
  assert result.lbfgs_error is None, result.lbfgs_error
  # no better point on a coarse grid around the result
  evaluator = llgi_scatfrac.llgi_scatfrac_b_factor_target_evaluator(
    f_eff=args["f_eff"], selection=args["working_selection"],
    f_calc=args["f_calc"], dobs=args["dobs"], sigmaa=args["sigmaa"],
    teps=args["teps"], resn=args["resn"], ss=args["d_star_sq"] / 4,
    centric_flags=args["centric_flags"], scale_factor=k,
    scatfrac_inf_start=result.scatfrac_inf, b_scatfrac_start=result.b_scatfrac,
    restraint_sigma=0.0, max_iterations=0,
    hybrid=llgi_hybrid.get_hybrid(llgi_data).with_rice_kappa(0.0))
  f_result, _ = evaluator.compute_functional_and_gradients()
  for dz in (-0.2, 0, 0.2):
    for db in (-5, 0, 5):
      evaluator.x = evaluator.x_from(result.scatfrac_inf * math.exp(dz),
        result.b_scatfrac + db)
      f, _ = evaluator.compute_functional_and_gradients()
      assert f >= f_result - 1.e-9, (dz, db, f, f_result)

def exercise_moments_form():
  """ scatfrac_moments recovers a resolution-dependent normalisation, and
  ScatFrac only rescales Ecalc in the F-scale functor: the target with
  ScatFrac S and scale k equals the target on k*Fcalc/sqrt(S) with
  ScatFrac 1 and k 1 (the variance depends on Dobs*sigmaA alone), with
  the gradient with respect to Fcalc by the chain rule. """
  from cctbx.xray import ext as xray_ext
  rs = random.Random(7)
  n = 4000
  dss = random_d_star_sq(n, rs)
  s_true = 0.6 + 2.0 * dss   # model normalisation rising with resolution
  resn = flex.double([rs.uniform(5, 15) for i in range(n)])
  f_calc = flex.complex_double()
  for i in range(n):
    sd = resn[i] * math.sqrt(s_true[i] / 2)
    f_calc.append(complex(rs.gauss(0, sd), rs.gauss(0, sd)))
  s_m = llgi_scatfrac.scatfrac_moments(f_calc, flex.double(n, 1.), resn, dss,
    n_bins=10)
  assert flex.max(flex.abs(s_m / s_true - 1)) < 0.25
  assert abs(flex.mean(s_m / s_true) - 1) < 0.03
  k = 1.3
  s_k = llgi_scatfrac.scatfrac_moments(f_calc, flex.double(n, 1.), resn, dss,
    n_bins=10, scale_factor=k)
  assert approx_equal(s_k, s_m * k * k)
  f_eff = flex.double([rs.uniform(1, 30) for i in range(n)])
  dobs = flex.double([rs.uniform(0.5, 1) for i in range(n)])
  sigmaa = flex.double(n, 0.8)
  common = dict(f_eff=f_eff, r_free_flags=flex.bool(n, False), dobs=dobs,
    sigmaa=sigmaa, teps=flex.double(n, 1.), resn=resn,
    centric_flags=flex.bool(n, False), compute_gradients=True)
  a = xray_ext.llgi_target_and_gradients(f_calc=f_calc, scatfrac=s_k,
    scale_factor=k, **common)
  factor = k / flex.sqrt(s_k)
  b = xray_ext.llgi_target_and_gradients(f_calc=f_calc * factor,
    scatfrac=flex.double(n, 1.), scale_factor=1., **common)
  assert approx_equal(a.target_work(), b.target_work(), eps=1e-12)
  assert flex.max(flex.abs(a.gradients_work() - b.gradients_work() * factor)) < 1e-12
  # the estimator returns the moments ScatFrac (the default method)
  params = llgi_scatfrac.llgi_scatfrac_params.extract()
  assert params.method == "moments"
  r = llgi_scatfrac.estimate_llgi_scatfrac(f_eff=f_eff,
    working_selection=flex.bool(n, True), f_calc=f_calc, dobs=dobs,
    sigmaa=sigmaa, teps=flex.double(n, 1.), resn=resn,
    centric_flags=flex.bool(n, False), d_star_sq=dss, scale_factor=k,
    params=params)
  params.moments_n_bins = 10
  r10 = llgi_scatfrac.estimate_llgi_scatfrac(f_eff=f_eff,
    working_selection=flex.bool(n, True), f_calc=f_calc, dobs=dobs,
    sigmaa=sigmaa, teps=flex.double(n, 1.), resn=resn,
    centric_flags=flex.bool(n, False), d_star_sq=dss, scale_factor=k,
    params=params)
  assert approx_equal(r10.scatfrac, s_k)
  assert r.scatfrac_inf is None

def exercise():
  exercise_moments_form()
  exercise_ratio_of_sums_robust_to_single_outlier_reflection()
  exercise_scatfrac_likelihood_uses_working_set_selection()
  exercise_scatfrac_b_factor_evaluator_gradient_finite_difference()
  exercise_scatfrac_b_factor_recovers_falling_trend()
  exercise_scatfrac_b_factor_recovers_rising_trend()
  exercise_b_factor_restraint_penalty_finite_difference()
  exercise_scatfrac_b_factor_restraint_pulls_toward_zero()
  exercise_likelihood_agrees_with_moments()
  exercise_scatfrac_fit_with_hybrid_converges()
  print("OK")

if (__name__ == "__main__"):
  exercise()
