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

def exercise():
  exercise_recovers_known_curves()
  exercise_convenience_wrapper_matches_two_step()
  exercise_scatfrac_uses_full_set_sigmaa_uses_test_set_only()
  print("OK")

if (__name__ == "__main__"):
  exercise()
