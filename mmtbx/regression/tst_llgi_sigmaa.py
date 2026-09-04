from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import mmtbx.refinement.llgi_sigmaa as llgi_sigmaa
from libtbx.test_utils import approx_equal
import math
import random

def rice_sample(ec, v, random_state):
  """ Draw one sample from the (acentric) Rice distribution with location
  ec (>= 0) and total variance v, matching the statistical model llgi.h's
  target_one_h assumes: the true complex structure factor sits at
  (ec, 0) (phase folded to zero WLOG) plus isotropic 2D Gaussian noise of
  total variance v (so each of Re/Im has variance v/2); the observed
  amplitude is the modulus of that noisy complex value. """
  sigma = math.sqrt(max(v, 1.e-10) / 2.0)
  re = random_state.gauss(ec, sigma)
  im = random_state.gauss(0.0, sigma)
  return math.sqrt(re * re + im * im)

def build_synthetic_dataset(n_refl, seed):
  """ Builds a self-consistent synthetic dataset for a known sigmaA(d*^2)
  and ScatFrac(d*^2) curve pair: Fcalc magnitudes are drawn from an
  exponential (Wilson-like acentric intensity) distribution with mean
  ScatFrac(d*^2)*Teps*Resn^2 (matching estimate_llgi_scatfrac's
  assumption exactly), and Feff is drawn from the true Rice distribution
  implied by D = Dobs*sigmaA/sqrt(ScatFrac), Ecalc = D*Fcalc, and
  V = Teps*Resn^2*(Teps - D^2) -- i.e. genuinely consistent with what
  cctbx.xray.targets.llgi's target_one_h assumes, not an approximation
  (an earlier draft of this test used abs(Normal(...)) as a stand-in for
  the Rice distribution; that biases the mean upward at low EC/sqrt(V)
  well beyond the Rice distribution's own bias, especially at high
  resolution where sigmaA is small, and was traced as the actual cause of
  an earlier, spurious-looking large recovery error before being
  replaced with a true Rice sampler here).
  """
  random_state = random.Random(seed)
  d_star_sq = flex.double(sorted(
    random_state.uniform(0.001, 0.25) for i in range(n_refl)))

  def true_sigmaa(dss):
    return 0.95 - 0.6 * dss / 0.25
  def true_scatfrac(dss):
    return 0.9 - 0.3 * dss / 0.25

  teps = flex.double(n_refl, 1.0)
  resn = flex.double(n_refl, 1.0)
  dobs = flex.double(n_refl, 0.85)
  sigmaa_true = flex.double([true_sigmaa(d) for d in d_star_sq])
  scatfrac_true = flex.double([true_scatfrac(d) for d in d_star_sq])

  f_calc_list = []
  f_eff_list = []
  for i in range(n_refl):
    mean_fc_sq = scatfrac_true[i] * teps[i] * resn[i] ** 2
    fc_mag = math.sqrt(random_state.expovariate(1.0 / mean_fc_sq))
    phase = random_state.uniform(0.0, 2.0 * math.pi)
    f_calc_list.append(
      complex(fc_mag * math.cos(phase), fc_mag * math.sin(phase)))
    d = dobs[i] * sigmaa_true[i] / math.sqrt(scatfrac_true[i])
    ec = d * fc_mag
    v_e = teps[i] - d * d
    v = teps[i] * resn[i] ** 2 * v_e
    f_eff_list.append(rice_sample(ec, v, random_state))

  return dict(
    f_eff=flex.double(f_eff_list),
    f_calc=flex.complex_double(f_calc_list),
    dobs=dobs, teps=teps, resn=resn, d_star_sq=d_star_sq,
    centric_flags=flex.bool(n_refl, False),
    sigmaa_true=sigmaa_true, scatfrac_true=scatfrac_true)

def exercise_recovers_known_curves():
  # With enough reflections and a physically self-consistent (true Rice-
  # distributed) synthetic dataset, the two-step estimator
  # (estimate_llgi_scatfrac then estimate_llgi_sigmaa) should recover
  # both curves reasonably closely -- this is the estimator's core
  # correctness property, and is a much stronger check than "does not
  # crash": an earlier, structurally different (jointly-optimised)
  # design was caught by exactly this kind of test failing to recover
  # the truth at all (oscillating, bound-hugging fits), traced to a real
  # non-identifiability between sigmaA and ScatFrac when both are fit
  # from the LLGI target alone -- see doc/llgi_target_design.md sec. 5.2.
  n_refl = 6000
  data = build_synthetic_dataset(n_refl, seed=7)
  r_free_flags = flex.bool([(i % 5 == 0) for i in range(n_refl)])

  scatfrac_fit = llgi_sigmaa.estimate_llgi_scatfrac(
    f_calc=data["f_calc"], teps=data["teps"], resn=data["resn"],
    d_star_sq=data["d_star_sq"], centric_flags=data["centric_flags"],
    n_coeffs=6)
  scatfrac_err = flex.mean(flex.abs(scatfrac_fit - data["scatfrac_true"]))
  assert scatfrac_err < 0.05, scatfrac_err

  result = llgi_sigmaa.estimate_llgi_sigmaa(
    f_eff=data["f_eff"], r_free_flags=r_free_flags, f_calc=data["f_calc"],
    dobs=data["dobs"], scatfrac=scatfrac_fit, teps=data["teps"],
    resn=data["resn"], centric_flags=data["centric_flags"],
    d_star_sq=data["d_star_sq"])
  sigmaa_err = flex.mean(flex.abs(result.sigmaa - data["sigmaa_true"]))
  # A loose bound: this is a statistical recovery test (finite-sample
  # LBFGS fit against noisy synthetic data), not an exact-answer check,
  # so the threshold is set well above typical run-to-run sampling
  # variance rather than tuned tightly to one seed. The point of this
  # test is to catch a fit that is badly wrong (e.g. the ~0.3-0.5 mean
  # abs errors, bound-hugging/oscillating fits seen during development
  # with a non-identifiable joint sigmaA/ScatFrac design, or with an
  # incorrectly-generated synthetic dataset -- see build_synthetic_
  # dataset's docstring), not to certify a specific numerical precision.
  assert sigmaa_err < 0.1, sigmaa_err

def exercise_convenience_wrapper_matches_two_step():
  # estimate_llgi_sigmaa_scatfrac should give identical results to
  # calling estimate_llgi_scatfrac then estimate_llgi_sigmaa by hand.
  n_refl = 500
  data = build_synthetic_dataset(n_refl, seed=11)
  r_free_flags = flex.bool([(i % 5 == 0) for i in range(n_refl)])

  scatfrac_fit = llgi_sigmaa.estimate_llgi_scatfrac(
    f_calc=data["f_calc"], teps=data["teps"], resn=data["resn"],
    d_star_sq=data["d_star_sq"], centric_flags=data["centric_flags"],
    n_coeffs=8)
  step2 = llgi_sigmaa.estimate_llgi_sigmaa(
    f_eff=data["f_eff"], r_free_flags=r_free_flags, f_calc=data["f_calc"],
    dobs=data["dobs"], scatfrac=scatfrac_fit, teps=data["teps"],
    resn=data["resn"], centric_flags=data["centric_flags"],
    d_star_sq=data["d_star_sq"])

  combined = llgi_sigmaa.estimate_llgi_sigmaa_scatfrac(
    f_eff=data["f_eff"], r_free_flags=r_free_flags, f_calc=data["f_calc"],
    dobs=data["dobs"], teps=data["teps"], resn=data["resn"],
    centric_flags=data["centric_flags"], d_star_sq=data["d_star_sq"])

  assert approx_equal(list(combined.scatfrac), list(scatfrac_fit))
  assert approx_equal(list(combined.sigmaa), list(step2.sigmaa))
  assert approx_equal(combined.target, step2.target)

def exercise_scatfrac_uses_full_set_sigmaa_uses_test_set_only():
  # estimate_llgi_scatfrac's result must not depend on r_free_flags at
  # all (it does not take r_free_flags as an argument -- this just
  # confirms the API shape matches the design intent: ScatFrac is fit on
  # the full reflection set, unlike sigmaA which is restricted to the
  # R-free/test set). A no-op check by construction, kept as a guard
  # against a future signature change accidentally reintroducing a
  # selection argument to estimate_llgi_scatfrac.
  import inspect
  argspec = inspect.getfullargspec(llgi_sigmaa.estimate_llgi_scatfrac)
  assert "r_free_flags" not in argspec.args
  argspec2 = inspect.getfullargspec(llgi_sigmaa.estimate_llgi_sigmaa)
  assert "r_free_flags" in argspec2.args

def exercise_scatfrac_robust_to_single_outlier_reflection():
  # Real bug, found running target=llgi against real (2g38) data: a
  # single reflection whose |Fcalc| happened to swing ~10x between two
  # ordinary macrocycles of refinement (small Resn amplified its squared
  # contribution) inflated an entire ~1600-reflection high-resolution
  # bin's ScatFrac estimate by ~10x on its own -- traced to
  # estimate_llgi_scatfrac averaging the per-reflection ratio
  # Fcalc_i^2/(Teps_i*Resn_i^2) directly (mean(X/Y)) rather than summing
  # numerator and denominator separately before dividing (ratio of
  # sums, equivalently ratio of means) -- mathematically different
  # quantities for a heavy-tailed numerator, and the mean-of-ratios
  # version is far more sensitive to a single extreme value.
  #
  # Getting this test to actually reproduce the failure took two tries:
  # a first attempt used a CONSTANT Resn across the bin, under which
  # mean(X/Y) and ratio-of-sums(X,Y) are algebraically identical (both
  # reduce to mean(X)/constant) -- so it could not distinguish the fixed
  # code from the buggy code at all, and passed even before the fix was
  # applied. The real bug specifically needs Resn to vary substantially
  # WITHIN one resolution bin (confirmed on real data: Resn ranged
  # ~1.97-183 within a single ~1600-reflection high-resolution bin, a
  # ~93x spread, std/mean ~0.86 -- resolution alone does not pin down
  # Resn tightly) and the outlier to sit at unusually SMALL Resn (as it
  # did in the real case) -- only then does averaging the ratio directly
  # (dividing by that reflection's own tiny Resn^2 before any averaging)
  # diverge from summing numerator and denominator separately.
  n_refl = 200
  random_state = random.Random(13)
  d_star_sq = flex.double([0.22 + 0.001 * i for i in range(n_refl)])
  teps = flex.double(n_refl, 1.0)
  # Resn varying log-uniformly across a ~90x range within the bin,
  # matching the real data's spread.
  resn_list = [
    math.exp(random_state.uniform(math.log(2.0), math.log(180.0)))
    for i in range(n_refl)]
  true_scatfrac = 0.6
  f_calc_list = []
  for i in range(n_refl):
    mean_fc_sq_i = true_scatfrac * teps[i] * resn_list[i] ** 2
    fc_mag = math.sqrt(random_state.expovariate(1.0 / mean_fc_sq_i))
    phase = random_state.uniform(0.0, 2.0 * math.pi)
    f_calc_list.append(complex(fc_mag * math.cos(phase), fc_mag * math.sin(phase)))
  # Implant one outlier at the SMALLEST Resn in the bin (matching the
  # real failure exactly), with a moderately (10x sqrt-scale) elevated
  # |Fcalc| relative to what that reflection's own small Resn would
  # typically produce.
  small_resn_idx = min(range(n_refl), key=lambda i: resn_list[i])
  mean_fc_sq_outlier = true_scatfrac * teps[small_resn_idx] * resn_list[small_resn_idx] ** 2
  f_calc_list[small_resn_idx] = complex(10.0 * math.sqrt(mean_fc_sq_outlier), 0.0)
  f_calc = flex.complex_double(f_calc_list)
  resn = flex.double(resn_list)
  centric_flags = flex.bool(n_refl, False)

  scatfrac_fit = llgi_sigmaa.estimate_llgi_scatfrac(
    f_calc=f_calc, teps=teps, resn=resn, d_star_sq=d_star_sq,
    centric_flags=centric_flags, n_coeffs=4, spline_degree=3, n_bins=1)
  fitted = flex.mean(scatfrac_fit)
  # A single outlier reflection (even a large one) at small Resn should
  # not be able to drag the ratio-of-sums fit far from the true value --
  # this specific construction reproduces a ~50% relative overestimate
  # under the pre-fix mean-of-ratios implementation (empirically
  # checked against a standalone reimplementation of both formulas
  # during development of this fix), well outside this bound.
  assert abs(fitted - true_scatfrac) < 0.15, (fitted, true_scatfrac)

def exercise_curvature_penalty_finite_difference():
  # Finite-difference check of _spline_curvature_penalty_and_gradient
  # itself (the analytic second-difference gradient), independent of any
  # LLGI machinery.
  import numpy as np
  random_state = np.random.RandomState(23)
  for trial in range(5):
    n = random_state.randint(3, 10)
    coeffs = random_state.normal(size=n)
    weight = random_state.uniform(0.01, 2.0)
    f0, g0 = llgi_sigmaa._spline_curvature_penalty_and_gradient(
      coeffs, weight)
    eps = 1.e-6
    g_fd = np.zeros_like(coeffs)
    for i in range(n):
      cp = coeffs.copy(); cp[i] += eps
      cm = coeffs.copy(); cm[i] -= eps
      fp, _ = llgi_sigmaa._spline_curvature_penalty_and_gradient(cp, weight)
      fm, _ = llgi_sigmaa._spline_curvature_penalty_and_gradient(cm, weight)
      g_fd[i] = (fp - fm) / (2 * eps)
    assert approx_equal(list(g0), list(g_fd), eps=1.e-6)

  # Edge cases: weight=0 and n<3 must be exact no-ops (used by callers
  # to disable the restraint entirely without special-casing).
  f_w0, g_w0 = llgi_sigmaa._spline_curvature_penalty_and_gradient(
    np.array([1.0, -2.0, 3.0, 0.5]), 0.0)
  assert f_w0 == 0.0 and list(g_w0) == [0.0, 0.0, 0.0, 0.0]
  f_n2, g_n2 = llgi_sigmaa._spline_curvature_penalty_and_gradient(
    np.array([1.0, -2.0]), 5.0)
  assert f_n2 == 0.0 and list(g_n2) == [0.0, 0.0]

  # An exactly-linear (or constant) coefficient vector has zero second
  # difference everywhere, so the penalty (and its gradient) must vanish
  # regardless of weight -- this is the property that makes the restraint
  # "invisible" on an already log-linear-like fit.
  c_line = np.linspace(-3.0, 4.0, 8)
  f_line, g_line = llgi_sigmaa._spline_curvature_penalty_and_gradient(
    c_line, 10.0)
  assert approx_equal(f_line, 0.0, eps=1.e-10)
  assert approx_equal(list(g_line), [0.0] * 8, eps=1.e-10)

def exercise_curvature_penalty_evaluator_gradient_finite_difference():
  # End-to-end finite-difference check of
  # llgi_sigmaa_target_evaluator.compute_functional_and_gradients with
  # the curvature restraint switched ON (curvature_weight > 0), against
  # a small synthetic dataset -- catches a sign error or an off-by-one
  # in how the penalty's gradient is added onto the LLGI target's own
  # spline-coefficient gradient (as opposed to
  # exercise_curvature_penalty_finite_difference above, which only
  # checks the penalty term in isolation).
  import numpy as np
  n_refl = 400
  data = build_synthetic_dataset(n_refl, seed=29)
  r_free_flags = flex.bool([(i % 5 == 0) for i in range(n_refl)])
  scatfrac_fit = llgi_sigmaa.estimate_llgi_scatfrac(
    f_calc=data["f_calc"], teps=data["teps"], resn=data["resn"],
    d_star_sq=data["d_star_sq"], centric_flags=data["centric_flags"],
    n_coeffs=6)
  sigmaa_design = llgi_sigmaa._b_spline_design_matrix(
    data["d_star_sq"].as_numpy_array(), 6, 3)
  evaluator = llgi_sigmaa.llgi_sigmaa_target_evaluator(
    f_eff=data["f_eff"], selection=r_free_flags, f_calc=data["f_calc"],
    dobs=data["dobs"], scatfrac=scatfrac_fit, teps=data["teps"],
    resn=data["resn"], centric_flags=data["centric_flags"],
    scale_factor=1.0, sigmaa_design=sigmaa_design, n_sigmaa_coeffs=6,
    max_iterations=0,  # Don't actually optimise -- just probe the
                        # gradient at the LBFGS starting point (x=0).
    curvature_weight=0.35)

  x0 = np.array(evaluator.x)
  f0, g0 = evaluator.compute_functional_and_gradients()
  g0 = np.array(g0)
  eps = 1.e-6
  g_fd = np.zeros_like(x0)
  for i in range(len(x0)):
    xp = x0.copy(); xp[i] += eps
    evaluator.x = flex.double(xp)
    fp, _ = evaluator.compute_functional_and_gradients()
    xm = x0.copy(); xm[i] -= eps
    evaluator.x = flex.double(xm)
    fm, _ = evaluator.compute_functional_and_gradients()
    g_fd[i] = (fp - fm) / (2 * eps)
  evaluator.x = flex.double(x0)
  assert approx_equal(list(g0), list(g_fd), eps=1.e-3)

def exercise_curvature_penalty_tames_synthetic_collapse():
  # The restraint's actual purpose: on a dataset where a sparse, noisy
  # high-resolution tail lets the LLGI likelihood alone drive sigmaA
  # sharply toward its lower bound in the last shell or two (mirroring
  # what was observed on real 2G38 data -- see the E-scale first-fit
  # diagnostic artifacts), a small curvature_weight should measurably
  # damp that collapse relative to curvature_weight=0, without
  # materially moving the well-determined (low/mid-resolution, many
  # reflections) part of the curve.
  #
  # Construct this directly: most of the resolution range has plenty of
  # reflections at a true sigmaA consistent with a smooth trend, but the
  # last shell has very few R-free reflections and its Feff values are
  # generated independently of any smooth trend (pure noise), so the
  # unrestrained LLGI fit is free to swing far from its neighbours
  # there.
  n_main = 1200
  n_tail = 15
  random_state = random.Random(41)

  d_star_sq_main = flex.double(sorted(
    random_state.uniform(0.001, 0.20) for i in range(n_main)))
  d_star_sq_tail = flex.double(sorted(
    random_state.uniform(0.205, 0.21) for i in range(n_tail)))
  d_star_sq = d_star_sq_main
  d_star_sq.extend(d_star_sq_tail)

  def true_sigmaa(dss):
    return 0.9 - 0.5 * dss / 0.21  # mild, smooth, log-plausible trend

  n_refl = n_main + n_tail
  teps = flex.double(n_refl, 1.0)
  resn = flex.double(n_refl, 1.0)
  dobs = flex.double(n_refl, 0.85)
  scatfrac_true = flex.double(n_refl, 0.75)

  f_calc_list = []
  f_eff_list = []
  for i in range(n_refl):
    mean_fc_sq = scatfrac_true[i] * teps[i] * resn[i] ** 2
    fc_mag = math.sqrt(random_state.expovariate(1.0 / mean_fc_sq))
    phase = random_state.uniform(0.0, 2.0 * math.pi)
    f_calc_list.append(
      complex(fc_mag * math.cos(phase), fc_mag * math.sin(phase)))
    if(i < n_main):
      sa = true_sigmaa(d_star_sq[i])
      d = dobs[i] * sa / math.sqrt(scatfrac_true[i])
      ec = d * fc_mag
      v_e = teps[i] - d * d
      v = teps[i] * resn[i] ** 2 * v_e
      f_eff_list.append(rice_sample(ec, v, random_state))
    else:
      # Tail reflections: Feff drawn independent of any smooth sigmaA
      # trend (just Wilson-like noise), so nothing in the data itself
      # anchors this shell to the main trend -- an unrestrained fit is
      # free to plunge toward the sigmoid floor here if doing so is even
      # marginally favoured by the (very sparse, noisy) test-set
      # likelihood in that shell.
      f_eff_list.append(math.sqrt(
        random_state.expovariate(1.0 / (teps[i] * resn[i] ** 2))))

  f_calc = flex.complex_double(f_calc_list)
  f_eff = flex.double(f_eff_list)
  centric_flags = flex.bool(n_refl, False)
  # Every tail reflection is R-free (test set), matching the "very few
  # informative reflections in the last shell" scenario this restraint
  # targets; keep the usual 1-in-5 split over the main range.
  r_free_flags = flex.bool(
    [(i % 5 == 0) for i in range(n_main)] + [True] * n_tail)

  common_kwargs = dict(
    f_eff=f_eff, r_free_flags=r_free_flags, f_calc=f_calc, dobs=dobs,
    scatfrac=scatfrac_true, teps=teps, resn=resn,
    centric_flags=centric_flags, d_star_sq=d_star_sq)

  params_unrestrained = llgi_sigmaa.llgi_sigmaa_scatfrac_params.extract()
  params_unrestrained.n_sigmaa_coeffs = 8
  params_unrestrained.sigmaa_curvature_weight = 0.0
  unrestrained = llgi_sigmaa.estimate_llgi_sigmaa(
    params=params_unrestrained, **common_kwargs)

  params_restrained = llgi_sigmaa.llgi_sigmaa_scatfrac_params.extract()
  params_restrained.n_sigmaa_coeffs = 8
  params_restrained.sigmaa_curvature_weight = 0.05
  restrained = llgi_sigmaa.estimate_llgi_sigmaa(
    params=params_restrained, **common_kwargs)

  # The last-shell sigmaA (evaluated at the tail reflections) should sit
  # closer to the smooth trend's extrapolated value with the restraint
  # on than with it off.
  true_tail = flex.double(
    [true_sigmaa(d_star_sq[i]) for i in range(n_main, n_refl)])
  err_unrestrained = flex.mean(
    flex.abs(unrestrained.sigmaa[n_main:] - true_tail))
  err_restrained = flex.mean(
    flex.abs(restrained.sigmaa[n_main:] - true_tail))
  assert err_restrained < err_unrestrained, (
    err_restrained, err_unrestrained)

  # And it should not have materially disturbed the LOW-resolution part
  # of the main range: this is a global 8-coefficient B-spline, so a
  # restraint strong enough to pull a genuinely severe collapse (the
  # unrestrained fit above hits the sigmoid floor, sigmaA=0.01, outright)
  # back toward the trend necessarily also reshapes some of the main
  # range near the tail, where the same basis functions that cover the
  # tail have non-negligible support -- that is expected, not a bug.
  # Far from the tail (the first fifth of the main range, d*^2 well
  # below the collapsing region), basis-function support for the tail is
  # negligible, so this restraint should leave that part of the curve
  # close to the unrestrained fit.
  main_diff_far = flex.mean(flex.abs(
    restrained.sigmaa[:n_main // 5] - unrestrained.sigmaa[:n_main // 5]))
  assert main_diff_far < 0.06, main_diff_far

def exercise_scatfrac_likelihood_evaluator_gradient_finite_difference():
  # End-to-end finite-difference check of llgi_scatfrac_target_evaluator.
  # compute_functional_and_gradients (the F-scale LLGI ScatFrac fit,
  # sigmaA held fixed) -- mirrors the equivalent sigmaA-evaluator check.
  # curvature_weight > 0 here specifically: this must exercise the
  # restraint's gradient mixed into the LLGI gradient, not just the
  # unrestrained path (which a curvature_weight=0.0 call would leave
  # entirely untested -- see exercise_curvature_penalty_finite_difference
  # for the restraint term checked in isolation).
  import numpy as np
  n_refl = 400
  data = build_synthetic_dataset(n_refl, seed=51)
  working_selection = flex.bool([(i % 5 != 0) for i in range(n_refl)])
  scatfrac_design = llgi_sigmaa._b_spline_design_matrix(
    data["d_star_sq"].as_numpy_array(), 6, 3)
  # A plausible starting curve (does not need to be the truth -- this is
  # only checking the gradient at whatever point LBFGS starts from).
  scatfrac_start = flex.double(n_refl, 0.6)
  evaluator = llgi_sigmaa.llgi_scatfrac_target_evaluator(
    f_eff=data["f_eff"], selection=working_selection, f_calc=data["f_calc"],
    dobs=data["dobs"], sigmaa=data["sigmaa_true"], teps=data["teps"],
    resn=data["resn"], centric_flags=data["centric_flags"],
    scale_factor=1.0, scatfrac_design=scatfrac_design,
    n_scatfrac_coeffs=6, scatfrac_start=scatfrac_start,
    max_iterations=0,  # probe the gradient at the LBFGS starting point
    curvature_weight=0.4)

  x0 = np.array(evaluator.x)
  f0, g0 = evaluator.compute_functional_and_gradients()
  g0 = np.array(g0)
  eps = 1.e-6
  g_fd = np.zeros_like(x0)
  for i in range(len(x0)):
    xp = x0.copy(); xp[i] += eps
    evaluator.x = flex.double(xp)
    fp, _ = evaluator.compute_functional_and_gradients()
    xm = x0.copy(); xm[i] -= eps
    evaluator.x = flex.double(xm)
    fm, _ = evaluator.compute_functional_and_gradients()
    g_fd[i] = (fp - fm) / (2 * eps)
  evaluator.x = flex.double(x0)
  assert approx_equal(list(g0), list(g_fd), eps=1.e-3)

def exercise_scatfrac_likelihood_recovers_known_curve_and_can_exceed_one():
  # estimate_llgi_scatfrac_likelihood, with sigmaA fixed at its TRUE
  # value, should recover the true ScatFrac(resolution) curve reasonably
  # well -- the E-then-F scheme's core correctness property, mirroring
  # exercise_recovers_known_curves for the historical F-then-... scheme.
  # Also checks the no-upper-bound-of-1 property directly: the synthetic
  # dataset here uses a ScatFrac curve that exceeds 1 (simulating Feff
  # not being on absolute scale -- see llgi_scatfrac_target_evaluator's
  # docstring), which estimate_llgi_scatfrac (the moment estimator) can
  # also represent, but which a sigmoid-bounded LBFGS reparameterisation
  # (like sigmaA's) could NOT.
  n_refl = 6000
  random_state = random.Random(53)
  d_star_sq = flex.double(sorted(
    random_state.uniform(0.001, 0.25) for i in range(n_refl)))

  def true_sigmaa(dss):
    return 0.7 - 0.3 * dss / 0.25
  def true_scatfrac(dss):
    # Deliberately exceeds 1 (up to ~1.3) at low resolution, simulating
    # an Feff absolute-scale offset (e.g. a wrong assumed solvent
    # fraction) rather than a "ScatFrac must stay below 1" ideal.
    return 1.3 - 0.7 * dss / 0.25

  teps = flex.double(n_refl, 1.0)
  resn = flex.double(n_refl, 1.0)
  dobs = flex.double(n_refl, 0.85)
  sigmaa_true = flex.double([true_sigmaa(d) for d in d_star_sq])
  scatfrac_true = flex.double([true_scatfrac(d) for d in d_star_sq])

  f_calc_list = []
  f_eff_list = []
  for i in range(n_refl):
    mean_fc_sq = scatfrac_true[i] * teps[i] * resn[i] ** 2
    fc_mag = math.sqrt(random_state.expovariate(1.0 / mean_fc_sq))
    phase = random_state.uniform(0.0, 2.0 * math.pi)
    f_calc_list.append(
      complex(fc_mag * math.cos(phase), fc_mag * math.sin(phase)))
    d = dobs[i] * sigmaa_true[i] / math.sqrt(scatfrac_true[i])
    ec = d * fc_mag
    v_e = teps[i] - d * d
    v = teps[i] * resn[i] ** 2 * v_e
    f_eff_list.append(rice_sample(ec, v, random_state))

  f_calc = flex.complex_double(f_calc_list)
  f_eff = flex.double(f_eff_list)
  centric_flags = flex.bool(n_refl, False)
  working_selection = flex.bool([(i % 5 != 0) for i in range(n_refl)])

  params = llgi_sigmaa.llgi_sigmaa_scatfrac_params.extract()
  params.n_scatfrac_coeffs = 6
  # The synthetic truth here is linear in d*^2 -- exactly representable by
  # the spline basis this test is actually exercising (see docstring:
  # "mirroring exercise_recovers_known_curves for the historical F-then-
  # ... scheme"), but NOT by the (now-default) b_factor mode, which is
  # log-linear in ss and would show a real, expected model-mismatch error
  # against this particular truth. Pin spline explicitly so this keeps
  # testing what it was written to test regardless of the phil default.
  params.scatfrac_model = "spline"
  result = llgi_sigmaa.estimate_llgi_scatfrac_likelihood(
    f_eff=f_eff, working_selection=working_selection, f_calc=f_calc,
    dobs=dobs, sigmaa=sigmaa_true, teps=teps, resn=resn,
    centric_flags=centric_flags, d_star_sq=d_star_sq, params=params)

  err = flex.mean(flex.abs(result.scatfrac - scatfrac_true))
  # Loose bound, matching exercise_recovers_known_curves's own rationale:
  # a statistical recovery test against noisy synthetic data, meant to
  # catch a badly wrong fit, not certify a specific precision.
  assert err < 0.15, err
  assert flex.max(result.scatfrac) > 1.0, flex.max(result.scatfrac)

def exercise_scatfrac_likelihood_uses_working_set_selection():
  # estimate_llgi_scatfrac_likelihood's result must actually depend on
  # working_selection (unlike estimate_llgi_scatfrac, the moment
  # estimator, which does not take a selection at all): fit the SAME
  # reflections but with two disjoint selections, each anchored on data
  # drawn from a DIFFERENT true ScatFrac curve, and confirm the two fits
  # land close to their own selection's truth -- not close to each other
  # or to some blend, which is what a selection-blind implementation
  # would produce.
  n_refl_each = 3000
  random_state = random.Random(56)

  def make_half(seed, scatfrac_level, d_lo=0.001, d_hi=0.25):
    d_star_sq = flex.double(sorted(
      random_state.uniform(d_lo, d_hi) for i in range(n_refl_each)))
    sigmaa_true = flex.double(n_refl_each, 0.6)
    scatfrac_true = flex.double(n_refl_each, scatfrac_level)
    teps = flex.double(n_refl_each, 1.0)
    resn = flex.double(n_refl_each, 1.0)
    dobs = flex.double(n_refl_each, 0.85)
    f_calc_list = []; f_eff_list = []
    for i in range(n_refl_each):
      mean_fc_sq = scatfrac_true[i] * teps[i] * resn[i] ** 2
      fc_mag = math.sqrt(random_state.expovariate(1.0 / mean_fc_sq))
      phase = random_state.uniform(0.0, 2.0 * math.pi)
      f_calc_list.append(
        complex(fc_mag * math.cos(phase), fc_mag * math.sin(phase)))
      d = dobs[i] * sigmaa_true[i] / math.sqrt(scatfrac_true[i])
      ec = d * fc_mag
      v = teps[i] * resn[i] ** 2 * (teps[i] - d * d)
      f_eff_list.append(rice_sample(ec, v, random_state))
    return dict(
      d_star_sq=d_star_sq, sigmaa_true=sigmaa_true,
      scatfrac_true=scatfrac_true, teps=teps, resn=resn, dobs=dobs,
      f_calc=flex.complex_double(f_calc_list),
      f_eff=flex.double(f_eff_list))

  low = make_half(56, scatfrac_level=0.4)
  high = make_half(56, scatfrac_level=1.1)  # deliberately > 1

  def cat(key):
    return low[key].concatenate(high[key])
  d_star_sq = cat("d_star_sq")
  f_calc = flex.complex_double(list(low["f_calc"]) + list(high["f_calc"]))
  f_eff = cat("f_eff")
  dobs = cat("dobs")
  sigmaa_true = cat("sigmaa_true")
  teps = cat("teps")
  resn = cat("resn")
  centric_flags = flex.bool(2 * n_refl_each, False)

  select_low = flex.bool(
    [True] * n_refl_each + [False] * n_refl_each)
  select_high = flex.bool(
    [False] * n_refl_each + [True] * n_refl_each)

  params = llgi_sigmaa.llgi_sigmaa_scatfrac_params.extract()
  params.n_scatfrac_coeffs = 4
  fit_low = llgi_sigmaa.estimate_llgi_scatfrac_likelihood(
    f_eff=f_eff, working_selection=select_low, f_calc=f_calc, dobs=dobs,
    sigmaa=sigmaa_true, teps=teps, resn=resn,
    centric_flags=centric_flags, d_star_sq=d_star_sq, params=params)
  fit_high = llgi_sigmaa.estimate_llgi_scatfrac_likelihood(
    f_eff=f_eff, working_selection=select_high, f_calc=f_calc, dobs=dobs,
    sigmaa=sigmaa_true, teps=teps, resn=resn,
    centric_flags=centric_flags, d_star_sq=d_star_sq, params=params)

  # Each fit, evaluated on ITS OWN selection's reflections, should sit
  # near that selection's own truth (0.4 or 1.1) -- not the other one,
  # and not their midpoint (~0.75).
  mean_low_on_low = flex.mean(fit_low.scatfrac.select(select_low))
  mean_high_on_high = flex.mean(fit_high.scatfrac.select(select_high))
  assert abs(mean_low_on_low - 0.4) < 0.15, mean_low_on_low
  assert abs(mean_high_on_high - 1.1) < 0.15, mean_high_on_high
  assert mean_high_on_high > mean_low_on_low + 0.3, (
    mean_low_on_low, mean_high_on_high)

  try:
    llgi_sigmaa.estimate_llgi_scatfrac_likelihood(
      f_eff=f_eff, working_selection=flex.bool(2 * n_refl_each, False),
      f_calc=f_calc, dobs=dobs, sigmaa=sigmaa_true, teps=teps, resn=resn,
      centric_flags=centric_flags, d_star_sq=d_star_sq)
  except RuntimeError as e:
    assert "working-set" in str(e)
  else:
    raise RuntimeError("Expected RuntimeError for an empty working set.")

def exercise_scatfrac_curvature_penalty_tames_synthetic_wobble():
  # The restraint's actual purpose (added after real 2G38 evidence -- see
  # llgi_scatfrac_target_evaluator's docstring): a narrow resolution band
  # with sparse, noisy working-set data can pull an UNrestrained ScatFrac
  # fit into a localised dip that the surrounding, well-determined curve
  # does not support. Unlike sigmaA (which collapses toward a sigmoid
  # floor), ScatFrac has no boundary to collapse toward, so the failure
  # mode here is a spurious LOCAL bump/dip in the interior of the curve,
  # not a one-sided runaway -- construct that directly: most of the
  # resolution range has plenty of working-set reflections consistent
  # with a smooth ScatFrac trend, but one narrow interior band has very
  # few, independent-of-the-trend (noise-only) reflections, giving the
  # unrestrained fit room to dip there.
  n_main = 1400
  n_wobble = 20
  random_state = random.Random(61)

  d_star_sq_main = flex.double(sorted(
    random_state.uniform(0.001, 0.10) for i in range(n_main // 2)))
  d_star_sq_wobble = flex.double(sorted(
    random_state.uniform(0.115, 0.125) for i in range(n_wobble)))
  d_star_sq_main2 = flex.double(sorted(
    random_state.uniform(0.14, 0.24) for i in range(n_main - n_main // 2)))
  d_star_sq = d_star_sq_main
  d_star_sq.extend(d_star_sq_wobble)
  d_star_sq.extend(d_star_sq_main2)

  def true_scatfrac(dss):
    return 0.9 - 0.5 * dss / 0.24  # mild, smooth trend, no wobble in truth

  n_refl = n_main + n_wobble
  teps = flex.double(n_refl, 1.0)
  resn = flex.double(n_refl, 1.0)
  dobs = flex.double(n_refl, 0.85)
  sigmaa_true = flex.double(n_refl, 0.6)  # fixed, as this evaluator assumes

  f_calc_list = []
  f_eff_list = []
  wobble_lo = 0.115
  wobble_hi = 0.125
  for i in range(n_refl):
    dss = d_star_sq[i]
    if(wobble_lo <= dss <= wobble_hi):
      # Wobble-band reflections: Fcalc/Feff drawn independent of the
      # smooth trend (pure Wilson-like noise), so nothing anchors this
      # narrow band to the surrounding curve.
      true_sf_i = random_state.uniform(0.3, 1.5)
    else:
      true_sf_i = true_scatfrac(dss)
    mean_fc_sq = true_sf_i * teps[i] * resn[i] ** 2
    fc_mag = math.sqrt(random_state.expovariate(1.0 / mean_fc_sq))
    phase = random_state.uniform(0.0, 2.0 * math.pi)
    f_calc_list.append(
      complex(fc_mag * math.cos(phase), fc_mag * math.sin(phase)))
    d = dobs[i] * sigmaa_true[i] / math.sqrt(true_sf_i)
    ec = d * fc_mag
    v_e = teps[i] - d * d
    v = teps[i] * resn[i] ** 2 * v_e
    f_eff_list.append(rice_sample(ec, v, random_state))

  f_calc = flex.complex_double(f_calc_list)
  f_eff = flex.double(f_eff_list)
  centric_flags = flex.bool(n_refl, False)
  working_selection = flex.bool(n_refl, True)  # all working set, none R-free

  common_kwargs = dict(
    f_eff=f_eff, working_selection=working_selection, f_calc=f_calc,
    dobs=dobs, sigmaa=sigmaa_true, teps=teps, resn=resn,
    centric_flags=centric_flags, d_star_sq=d_star_sq)

  # This test is specifically about the SPLINE mode's curvature restraint
  # (scatfrac_curvature_weight is ignored in scalar/b_factor -- see
  # llgi_sigmaa_scatfrac_params' help) -- pin it explicitly so it keeps
  # testing that regardless of the phil default.
  params_unrestrained = llgi_sigmaa.llgi_sigmaa_scatfrac_params.extract()
  params_unrestrained.scatfrac_model = "spline"
  params_unrestrained.n_scatfrac_coeffs = 8
  params_unrestrained.scatfrac_curvature_weight = 0.0
  unrestrained = llgi_sigmaa.estimate_llgi_scatfrac_likelihood(
    params=params_unrestrained, **common_kwargs)

  params_restrained = llgi_sigmaa.llgi_sigmaa_scatfrac_params.extract()
  params_restrained.scatfrac_model = "spline"
  params_restrained.n_scatfrac_coeffs = 8
  params_restrained.scatfrac_curvature_weight = 0.1
  restrained = llgi_sigmaa.estimate_llgi_scatfrac_likelihood(
    params=params_restrained, **common_kwargs)

  # The wobble band's fitted ScatFrac should sit closer to the smooth
  # trend's interpolated value with the restraint on than off.
  wobble_indices = [i for i in range(n_refl)
                     if wobble_lo <= d_star_sq[i] <= wobble_hi]
  true_wobble = flex.double(
    [true_scatfrac(d_star_sq[i]) for i in wobble_indices])
  unrestrained_wobble = flex.double(
    [unrestrained.scatfrac[i] for i in wobble_indices])
  restrained_wobble = flex.double(
    [restrained.scatfrac[i] for i in wobble_indices])
  err_unrestrained = flex.mean(flex.abs(unrestrained_wobble - true_wobble))
  err_restrained = flex.mean(flex.abs(restrained_wobble - true_wobble))
  assert err_restrained < err_unrestrained, (
    err_restrained, err_unrestrained)

  # And it should not have materially disturbed the well-determined main
  # range on the FAR side of the wobble band (the last fifth of the
  # d*^2 range, i.e. n_refl's tail, well past d_star_sq_wobble): a
  # global 8-coefficient spline necessarily reshapes some of the curve
  # near the perturbed region (confirmed directly: the region BETWEEN
  # the wobble band and this far tail does shift, same as the sigmaA
  # collapse test's own finding for its adjacent region), but should
  # leave the distant, opposite-side part close to the unrestrained fit.
  far_indices = list(range(n_refl))[-(n_refl // 5):]
  main_diff_far = flex.mean(flex.abs(
    flex.double([restrained.scatfrac[i] for i in far_indices]) -
    flex.double([unrestrained.scatfrac[i] for i in far_indices])))
  assert main_diff_far < 0.1, main_diff_far

def exercise_scatfrac_scalar_evaluator_gradient_finite_difference():
  # End-to-end finite-difference check of
  # llgi_scatfrac_scalar_target_evaluator.compute_functional_and_gradients
  # (a single shared ScatFrac value, not a spline).
  import numpy as np
  n_refl = 400
  data = build_synthetic_dataset(n_refl, seed=71)
  working_selection = flex.bool([(i % 5 != 0) for i in range(n_refl)])
  evaluator = llgi_sigmaa.llgi_scatfrac_scalar_target_evaluator(
    f_eff=data["f_eff"], selection=working_selection, f_calc=data["f_calc"],
    dobs=data["dobs"], sigmaa=data["sigmaa_true"], teps=data["teps"],
    resn=data["resn"], centric_flags=data["centric_flags"],
    scale_factor=1.0, scatfrac_start=0.6,
    max_iterations=0)  # probe the gradient at the LBFGS starting point

  x0 = np.array(evaluator.x)
  f0, g0 = evaluator.compute_functional_and_gradients()
  g0 = np.array(g0)
  eps = 1.e-6
  xp = x0.copy(); xp[0] += eps
  evaluator.x = flex.double(xp)
  fp, _ = evaluator.compute_functional_and_gradients()
  xm = x0.copy(); xm[0] -= eps
  evaluator.x = flex.double(xm)
  fm, _ = evaluator.compute_functional_and_gradients()
  g_fd = (fp - fm) / (2 * eps)
  evaluator.x = flex.double(x0)
  assert approx_equal(g0[0], g_fd, eps=1.e-3)

def exercise_scatfrac_scalar_recovers_known_constant_and_can_exceed_one():
  # estimate_llgi_scatfrac_likelihood(params.scatfrac_model="scalar"), with
  # sigmaA fixed at its true value, should recover a known CONSTANT
  # ScatFrac (no resolution dependence in the synthetic truth) closely,
  # and the returned array should be genuinely constant across
  # reflections (not merely close to it) since there is only one
  # parameter. Also checks the no-upper-bound-of-1 property, same as the
  # spline version's equivalent test.
  n_refl = 4000
  random_state = random.Random(73)
  d_star_sq = flex.double(sorted(
    random_state.uniform(0.001, 0.25) for i in range(n_refl)))
  sigmaa_true = flex.double(n_refl, 0.6)
  true_scatfrac = 1.15  # deliberately > 1
  teps = flex.double(n_refl, 1.0)
  resn = flex.double(n_refl, 1.0)
  dobs = flex.double(n_refl, 0.85)

  f_calc_list = []
  f_eff_list = []
  for i in range(n_refl):
    mean_fc_sq = true_scatfrac * teps[i] * resn[i] ** 2
    fc_mag = math.sqrt(random_state.expovariate(1.0 / mean_fc_sq))
    phase = random_state.uniform(0.0, 2.0 * math.pi)
    f_calc_list.append(
      complex(fc_mag * math.cos(phase), fc_mag * math.sin(phase)))
    d = dobs[i] * sigmaa_true[i] / math.sqrt(true_scatfrac)
    ec = d * fc_mag
    v = teps[i] * resn[i] ** 2 * (teps[i] - d * d)
    f_eff_list.append(rice_sample(ec, v, random_state))

  f_calc = flex.complex_double(f_calc_list)
  f_eff = flex.double(f_eff_list)
  centric_flags = flex.bool(n_refl, False)
  working_selection = flex.bool(n_refl, True)

  params = llgi_sigmaa.llgi_sigmaa_scatfrac_params.extract()
  params.scatfrac_model = "scalar"
  result = llgi_sigmaa.estimate_llgi_scatfrac_likelihood(
    f_eff=f_eff, working_selection=working_selection, f_calc=f_calc,
    dobs=dobs, sigmaa=sigmaa_true, teps=teps, resn=resn,
    centric_flags=centric_flags, d_star_sq=d_star_sq, params=params)

  # Genuinely constant (single parameter -> identical value everywhere).
  assert flex.max(result.scatfrac) == flex.min(result.scatfrac)
  fitted = result.scatfrac[0]
  assert abs(fitted - true_scatfrac) < 0.1, (fitted, true_scatfrac)
  assert fitted > 1.0, fitted

def exercise_scatfrac_scalar_ignores_resolution_dependence():
  # With a synthetic truth that DOES vary strongly with resolution, the
  # scalar fit should land somewhere between the low- and high-resolution
  # extremes (a single compromise value), NOT track the trend -- this is
  # the scalar mode's known, accepted limitation (see llgi_scatfrac_
  # scalar_target_evaluator's docstring), verified directly so a future
  # change that accidentally made this mode resolution-sensitive again
  # would be caught.
  n_refl = 4000
  random_state = random.Random(77)
  d_star_sq = flex.double(sorted(
    random_state.uniform(0.001, 0.25) for i in range(n_refl)))
  sigmaa_true = flex.double(n_refl, 0.6)

  def true_scatfrac(dss):
    return 1.4 - 1.0 * dss / 0.25  # strong trend: 1.4 -> 0.4

  teps = flex.double(n_refl, 1.0)
  resn = flex.double(n_refl, 1.0)
  dobs = flex.double(n_refl, 0.85)

  f_calc_list = []
  f_eff_list = []
  for i in range(n_refl):
    sf_i = true_scatfrac(d_star_sq[i])
    mean_fc_sq = sf_i * teps[i] * resn[i] ** 2
    fc_mag = math.sqrt(random_state.expovariate(1.0 / mean_fc_sq))
    phase = random_state.uniform(0.0, 2.0 * math.pi)
    f_calc_list.append(
      complex(fc_mag * math.cos(phase), fc_mag * math.sin(phase)))
    d = dobs[i] * sigmaa_true[i] / math.sqrt(sf_i)
    ec = d * fc_mag
    v = teps[i] * resn[i] ** 2 * (teps[i] - d * d)
    f_eff_list.append(rice_sample(ec, v, random_state))

  f_calc = flex.complex_double(f_calc_list)
  f_eff = flex.double(f_eff_list)
  centric_flags = flex.bool(n_refl, False)
  working_selection = flex.bool(n_refl, True)

  params = llgi_sigmaa.llgi_sigmaa_scatfrac_params.extract()
  params.scatfrac_model = "scalar"
  result = llgi_sigmaa.estimate_llgi_scatfrac_likelihood(
    f_eff=f_eff, working_selection=working_selection, f_calc=f_calc,
    dobs=dobs, sigmaa=sigmaa_true, teps=teps, resn=resn,
    centric_flags=centric_flags, d_star_sq=d_star_sq, params=params)

  fitted = result.scatfrac[0]
  # Strictly between the two extremes, not equal to either -- a single
  # compromise value, as expected for this mode's known limitation.
  assert 0.4 < fitted < 1.4, fitted

def _make_ss(d_star_sq):
  # d_star_sq is always a flex.double at every call site in this file.
  return d_star_sq.as_numpy_array() / 4.0

def exercise_scatfrac_b_factor_evaluator_gradient_finite_difference():
  # End-to-end finite-difference check of llgi_scatfrac_b_factor_target_
  # evaluator.compute_functional_and_gradients (the two-parameter
  # ScatFrac_inf*exp(-B_scatfrac*ss) fit).
  import numpy as np
  n_refl = 400
  data = build_synthetic_dataset(n_refl, seed=81)
  working_selection = flex.bool([(i % 5 != 0) for i in range(n_refl)])
  ss = _make_ss(data["d_star_sq"])
  evaluator = llgi_sigmaa.llgi_scatfrac_b_factor_target_evaluator(
    f_eff=data["f_eff"], selection=working_selection, f_calc=data["f_calc"],
    dobs=data["dobs"], sigmaa=data["sigmaa_true"], teps=data["teps"],
    resn=data["resn"], ss=ss, centric_flags=data["centric_flags"],
    scale_factor=1.0, scatfrac_inf_start=0.7, b_scatfrac_start=15.0,
    max_iterations=0)  # probe the gradient at the LBFGS starting point

  x0 = np.array(evaluator.x)
  f0, g0 = evaluator.compute_functional_and_gradients()
  g0 = np.array(g0)
  eps = 1.e-6
  g_fd = np.zeros_like(x0)
  for i in range(len(x0)):
    xp = x0.copy(); xp[i] += eps
    evaluator.x = flex.double(xp)
    fp, _ = evaluator.compute_functional_and_gradients()
    xm = x0.copy(); xm[i] -= eps
    evaluator.x = flex.double(xm)
    fm, _ = evaluator.compute_functional_and_gradients()
    g_fd[i] = (fp - fm) / (2 * eps)
  evaluator.x = flex.double(x0)
  assert approx_equal(list(g0), list(g_fd), eps=1.e-3)

def _synthetic_b_factor_scatfrac_dataset(
      n_refl, seed, scatfrac_inf, b_scatfrac, d_lo=0.001, d_hi=0.25):
  # Shared helper: build a self-consistent synthetic dataset with a KNOWN
  # ScatFrac_inf*exp(-B_scatfrac*ss) truth (positive OR negative
  # B_scatfrac), fixed sigmaA, Rice-sampled Feff -- same construction
  # pattern as build_synthetic_dataset, specialised to a log-linear
  # (rather than flat or spline) ScatFrac truth.
  random_state = random.Random(seed)
  d_star_sq = flex.double(sorted(
    random_state.uniform(d_lo, d_hi) for i in range(n_refl)))
  sigmaa_true = flex.double(n_refl, 0.6)
  teps = flex.double(n_refl, 1.0)
  resn = flex.double(n_refl, 1.0)
  dobs = flex.double(n_refl, 0.85)

  def true_scatfrac(dss):
    ss = dss / 4.0
    return scatfrac_inf * math.exp(-b_scatfrac * ss)

  f_calc_list = []
  f_eff_list = []
  for i in range(n_refl):
    sf_i = true_scatfrac(d_star_sq[i])
    mean_fc_sq = sf_i * teps[i] * resn[i] ** 2
    fc_mag = math.sqrt(random_state.expovariate(1.0 / mean_fc_sq))
    phase = random_state.uniform(0.0, 2.0 * math.pi)
    f_calc_list.append(
      complex(fc_mag * math.cos(phase), fc_mag * math.sin(phase)))
    d = dobs[i] * sigmaa_true[i] / math.sqrt(sf_i)
    ec = d * fc_mag
    v = teps[i] * resn[i] ** 2 * (teps[i] - d * d)
    f_eff_list.append(rice_sample(ec, v, random_state))

  return dict(
    d_star_sq=d_star_sq, sigmaa_true=sigmaa_true, teps=teps, resn=resn,
    dobs=dobs, f_calc=flex.complex_double(f_calc_list),
    f_eff=flex.double(f_eff_list),
    centric_flags=flex.bool(n_refl, False),
    working_selection=flex.bool(n_refl, True))

def _assert_b_factor_curve_matches_truth_over_observed_range(
      result, d_star_sq, scatfrac_inf_true, b_scatfrac_true):
  # ScatFrac_inf is the z_inf -> d*^2 = 0 EXTRAPOLATED intercept -- with
  # synthetic data confined to a finite resolution range (as any real
  # dataset is), (ScatFrac_inf, B_scatfrac) trade off against each other
  # much like the classic Wilson-plot scale/B-factor ambiguity: many
  # (ScatFrac_inf, B_scatfrac) pairs give nearly the SAME curve over the
  # observed range while extrapolating to very different d*^2=0 values.
  # Caught directly during development: a fit recovered scatfrac_inf=0.50
  # against a true 0.95 (b_scatfrac 10.5 against a true 25.0) YET matched
  # the true curve to within 0.13 mean absolute / 22% mean relative
  # difference over the actual observed d*^2 range, with log-log
  # correlation 0.99999997 -- i.e. the right curve SHAPE, offset by
  # exactly the kind of extrapolation ambiguity a finite resolution range
  # cannot resolve. So this checks the fitted curve against the truth
  # over the OBSERVED range (what is actually identifiable), not the
  # individual (ScatFrac_inf, B_scatfrac) values (which are not, in
  # general, tightly determined by a finite d*^2 range) -- and confirms
  # b_scatfrac's SIGN is still recovered correctly, since sign is a much
  # coarser property than magnitude and not sensitive to this ambiguity.
  import math
  import numpy as np
  n = d_star_sq.size()
  true_vals = np.array([
    scatfrac_inf_true * math.exp(-b_scatfrac_true * (d_star_sq[i] / 4.0))
    for i in range(n)])
  fit_vals = np.array([result.scatfrac[i] for i in range(n)])
  mean_rel_diff = np.mean(np.abs(true_vals - fit_vals) / true_vals)
  assert mean_rel_diff < 0.3, mean_rel_diff
  log_corr = np.corrcoef(np.log(true_vals), np.log(fit_vals))[0, 1]
  assert log_corr > 0.999, log_corr
  same_sign = (result.b_scatfrac > 0) == (b_scatfrac_true > 0)
  assert same_sign, (result.b_scatfrac, b_scatfrac_true)

def exercise_scatfrac_b_factor_recovers_falling_trend():
  # Ordinary case: ScatFrac falls with resolution (positive B_scatfrac,
  # like a normal crystallographic B-factor), sigmaA fixed at truth.
  data = _synthetic_b_factor_scatfrac_dataset(
    n_refl=6000, seed=83, scatfrac_inf=0.95, b_scatfrac=25.0)
  params = llgi_sigmaa.llgi_sigmaa_scatfrac_params.extract()
  params.scatfrac_model = "b_factor"
  # Truth's |B_scatfrac|=25 is itself far from 0 -- pin the restraint off
  # so this checks the unrestrained evaluator's recovery of a genuine
  # trend, not the (default-on) restraint's deliberate pull toward 0.
  # See exercise_scatfrac_b_factor_restraint_pulls_toward_zero below for
  # the restraint's own behaviour.
  params.scatfrac_b_factor_restraint_sigma = 0.0
  result = llgi_sigmaa.estimate_llgi_scatfrac_likelihood(
    f_eff=data["f_eff"], working_selection=data["working_selection"],
    f_calc=data["f_calc"], dobs=data["dobs"], sigmaa=data["sigmaa_true"],
    teps=data["teps"], resn=data["resn"],
    centric_flags=data["centric_flags"], d_star_sq=data["d_star_sq"],
    params=params)
  assert result.scatfrac_inf is not None and result.b_scatfrac is not None
  _assert_b_factor_curve_matches_truth_over_observed_range(
    result, data["d_star_sq"], 0.95, 25.0)

def exercise_scatfrac_b_factor_recovers_rising_trend():
  # The user's motivating case: a partial model built from its best-
  # ordered components first can have ScatFrac RISE toward high
  # resolution -- a NEGATIVE B_scatfrac. Confirms the sign is genuinely
  # unconstrained and gets recovered correctly, not just that the
  # falling case (the more common/expected direction) works.
  data = _synthetic_b_factor_scatfrac_dataset(
    n_refl=6000, seed=89, scatfrac_inf=0.5, b_scatfrac=-25.0)
  params = llgi_sigmaa.llgi_sigmaa_scatfrac_params.extract()
  params.scatfrac_model = "b_factor"
  # See exercise_scatfrac_b_factor_recovers_falling_trend's comment: pin
  # the restraint off so this isolates the unrestrained evaluator's
  # ability to recover a genuine (here rising) trend.
  params.scatfrac_b_factor_restraint_sigma = 0.0
  result = llgi_sigmaa.estimate_llgi_scatfrac_likelihood(
    f_eff=data["f_eff"], working_selection=data["working_selection"],
    f_calc=data["f_calc"], dobs=data["dobs"], sigmaa=data["sigmaa_true"],
    teps=data["teps"], resn=data["resn"],
    centric_flags=data["centric_flags"], d_star_sq=data["d_star_sq"],
    params=params)
  _assert_b_factor_curve_matches_truth_over_observed_range(
    result, data["d_star_sq"], 0.5, -25.0)

  # And the actual fitted curve should rise with resolution (decrease
  # with d, i.e. increase with d*^2): check the fitted ScatFrac at the
  # highest- and lowest-resolution ends of the input.
  ds2 = data["d_star_sq"]
  order = sorted(range(ds2.size()), key=lambda i: ds2[i])
  low_res_val = result.scatfrac[order[0]]
  high_res_val = result.scatfrac[order[-1]]
  assert high_res_val > low_res_val, (low_res_val, high_res_val)

def exercise_b_factor_restraint_penalty_finite_difference():
  # _b_factor_restraint_penalty_and_gradient's own gradient, checked
  # against a plain finite difference, independent of the evaluator.
  sigma = 10.0
  for b in (-30.0, -1.0, 0.0, 4.0, 25.0):
    f0, g0 = llgi_sigmaa._b_factor_restraint_penalty_and_gradient(b, sigma)
    eps = 1.e-6
    fp, _ = llgi_sigmaa._b_factor_restraint_penalty_and_gradient(
      b + eps, sigma)
    fm, _ = llgi_sigmaa._b_factor_restraint_penalty_and_gradient(
      b - eps, sigma)
    g_fd = (fp - fm) / (2 * eps)
    assert approx_equal(g0, g_fd, eps=1.e-4)
  # sigma <= 0 disables the restraint (matches weight=0 elsewhere).
  assert llgi_sigmaa._b_factor_restraint_penalty_and_gradient(25.0, 0.0) \
    == (0.0, 0.0)

def exercise_scatfrac_b_factor_restraint_pulls_toward_zero():
  # Motivating case: sigmaA and ScatFrac trade off against each other
  # through the F-scale LLGI target (see llgi_scatfrac_b_factor_target_
  # evaluator's docstring), so a dataset that is well described by a
  # SMALL true B_scatfrac can still let an unrestrained fit wander to a
  # much larger |B_scatfrac| with little likelihood cost. Confirm the
  # default restraint (scatfrac_b_factor_restraint_sigma=10, matching
  # the phil default) pulls such a fit closer to 0 than leaving it off,
  # without flipping recovery of a genuine, comfortably-sized trend.
  data = _synthetic_b_factor_scatfrac_dataset(
    n_refl=6000, seed=97, scatfrac_inf=0.8, b_scatfrac=3.0)
  params_unrestrained = llgi_sigmaa.llgi_sigmaa_scatfrac_params.extract()
  params_unrestrained.scatfrac_model = "b_factor"
  params_unrestrained.scatfrac_b_factor_restraint_sigma = 0.0
  result_unrestrained = llgi_sigmaa.estimate_llgi_scatfrac_likelihood(
    f_eff=data["f_eff"], working_selection=data["working_selection"],
    f_calc=data["f_calc"], dobs=data["dobs"], sigmaa=data["sigmaa_true"],
    teps=data["teps"], resn=data["resn"],
    centric_flags=data["centric_flags"], d_star_sq=data["d_star_sq"],
    params=params_unrestrained)

  params_restrained = llgi_sigmaa.llgi_sigmaa_scatfrac_params.extract()
  params_restrained.scatfrac_model = "b_factor"
  assert params_restrained.scatfrac_b_factor_restraint_sigma == 10.0
  result_restrained = llgi_sigmaa.estimate_llgi_scatfrac_likelihood(
    f_eff=data["f_eff"], working_selection=data["working_selection"],
    f_calc=data["f_calc"], dobs=data["dobs"], sigmaa=data["sigmaa_true"],
    teps=data["teps"], resn=data["resn"],
    centric_flags=data["centric_flags"], d_star_sq=data["d_star_sq"],
    params=params_restrained)

  # The restraint should not have flipped the sign or badly damaged
  # recovery of a genuine trend this close to the restraint width.
  assert (result_restrained.b_scatfrac > 0) == \
    (result_unrestrained.b_scatfrac > 0)
  assert abs(result_restrained.b_scatfrac) < abs(
    result_unrestrained.b_scatfrac) + 1.e-6, (
      result_restrained.b_scatfrac, result_unrestrained.b_scatfrac)

def exercise():
  exercise_recovers_known_curves()
  exercise_convenience_wrapper_matches_two_step()
  exercise_scatfrac_uses_full_set_sigmaa_uses_test_set_only()
  exercise_scatfrac_robust_to_single_outlier_reflection()
  exercise_curvature_penalty_finite_difference()
  exercise_curvature_penalty_evaluator_gradient_finite_difference()
  exercise_curvature_penalty_tames_synthetic_collapse()
  exercise_scatfrac_likelihood_evaluator_gradient_finite_difference()
  exercise_scatfrac_likelihood_recovers_known_curve_and_can_exceed_one()
  exercise_scatfrac_likelihood_uses_working_set_selection()
  exercise_scatfrac_curvature_penalty_tames_synthetic_wobble()
  exercise_scatfrac_scalar_evaluator_gradient_finite_difference()
  exercise_scatfrac_scalar_recovers_known_constant_and_can_exceed_one()
  exercise_scatfrac_scalar_ignores_resolution_dependence()
  exercise_scatfrac_b_factor_evaluator_gradient_finite_difference()
  exercise_scatfrac_b_factor_recovers_falling_trend()
  exercise_scatfrac_b_factor_recovers_rising_trend()
  exercise_b_factor_restraint_penalty_finite_difference()
  exercise_scatfrac_b_factor_restraint_pulls_toward_zero()
  print("OK")

if (__name__ == "__main__"):
  exercise()
