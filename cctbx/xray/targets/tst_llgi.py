from __future__ import absolute_import, division, print_function
from cctbx.xray import ext
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from six.moves import range
import math

def make_inputs(n_refl, centric, seed=0):
  # Deterministic pseudo-random-ish inputs (no RNG dependency) covering a
  # spread of magnitudes, chosen so that dobs*sigmaa/sqrt(scatfrac) stays
  # comfortably below sqrt(teps) (i.e. v_e = teps - d^2 > 0) for every
  # reflection, avoiding the negative-variance guard in the basic checks.
  f_eff = flex.double([0.3 + 0.17 * ((seed + i) % 7) for i in range(n_refl)])
  dobs = flex.double([0.4 + 0.05 * ((seed + i) % 6) for i in range(n_refl)])
  sigmaa = flex.double([0.5 + 0.03 * ((seed + i) % 5) for i in range(n_refl)])
  scatfrac = flex.double([0.8 + 0.02 * ((seed + i) % 4) for i in range(n_refl)])
  teps = flex.double([1.0 for i in range(n_refl)])
  resn = flex.double([1.0 + 0.1 * ((seed + i) % 3) for i in range(n_refl)])
  centric_flags = flex.bool([centric for i in range(n_refl)])
  r_free_flags = flex.bool([(i % 4 == 0) for i in range(n_refl)])
  f_calc = flex.complex_double([
    complex(
      (0.6 + 0.11 * ((seed + i) % 5)) * math.cos(0.37 * i + 0.1 * seed),
      (0.6 + 0.11 * ((seed + i) % 5)) * math.sin(0.37 * i + 0.1 * seed))
    for i in range(n_refl)])
  return dict(
    f_eff=f_eff, r_free_flags=r_free_flags, f_calc=f_calc, dobs=dobs,
    sigmaa=sigmaa, scatfrac=scatfrac, teps=teps, resn=resn,
    centric_flags=centric_flags)

def target_work(inputs, scale_factor=1.0, compute_gradients=False):
  return ext.llgi_target_and_gradients(
    f_eff=inputs["f_eff"],
    r_free_flags=inputs["r_free_flags"],
    f_calc=inputs["f_calc"],
    dobs=inputs["dobs"],
    sigmaa=inputs["sigmaa"],
    scatfrac=inputs["scatfrac"],
    scale_factor=scale_factor,
    teps=inputs["teps"],
    resn=inputs["resn"],
    centric_flags=inputs["centric_flags"],
    compute_gradients=compute_gradients)

def exercise_finite_difference_gradients(centric):
  n_refl = 12
  inputs = make_inputs(n_refl, centric=centric, seed=3)
  result = target_work(inputs, compute_gradients=True)
  ana_grads = result.gradients_work()
  # gradients_work() only covers work reflections (r_free_flags == False);
  # map back to full-array indices to perturb the matching f_calc entry.
  work_indices = [i for i in range(n_refl) if not inputs["r_free_flags"][i]]
  assert len(work_indices) == len(ana_grads)
  eps = 1.e-6
  for k, ih in enumerate(work_indices):
    fin_grads = []
    for part in ("re", "im"):
      vals = []
      for signed_eps in (eps, -eps):
        f_calc_pert = inputs["f_calc"].deep_copy()
        c0 = f_calc_pert[ih]
        if part == "re":
          f_calc_pert[ih] = complex(c0.real + signed_eps, c0.imag)
        else:
          f_calc_pert[ih] = complex(c0.real, c0.imag + signed_eps)
        inputs_pert = dict(inputs)
        inputs_pert["f_calc"] = f_calc_pert
        t_pert = target_work(inputs_pert, compute_gradients=False)
        vals.append(t_pert.target_work())
      fin_grads.append((vals[0] - vals[1]) / (2 * eps))
    # gradients_work() values are already divided by n_work (see
    # target_and_gradients in llgi.h), matching target_work() which is
    # itself target_sum/n_work -- so no additional rescaling is needed.
    # The double conj() applied in llgi.h (once in d_target_one_h_over_fc
    # via conj(fc_complex), once again in target_and_gradients'
    # std::conj(...) wrapper) cancels, so gradients_work()'s real and
    # imaginary parts equal d(target_work)/d(Re fc) and d(target_work)/
    # d(Im fc) directly, confirmed empirically below.
    g = ana_grads[k]
    d_by_dre_ana = g.real
    d_by_dim_ana = g.imag
    assert approx_equal(d_by_dre_ana, fin_grads[0], eps=5.e-5), (
      centric, ih, d_by_dre_ana, fin_grads[0])
    assert approx_equal(d_by_dim_ana, fin_grads[1], eps=5.e-5), (
      centric, ih, d_by_dim_ana, fin_grads[1])

def exercise_resn_scale_invariance():
  # Scaling f_eff, f_calc and resn all by the same positive factor should
  # leave the target value exactly unchanged (RESN is a pure change-of-
  # variables normalization, see llgi.h docstring).
  inputs = make_inputs(8, centric=False, seed=1)
  base = target_work(inputs).target_work()
  scaled = dict(inputs)
  factor = 2.3
  scaled["f_eff"] = inputs["f_eff"] * factor
  scaled["f_calc"] = inputs["f_calc"] * factor
  scaled["resn"] = inputs["resn"] * factor
  scaled_val = target_work(scaled).target_work()
  assert approx_equal(base, scaled_val, eps=1.e-9)

def exercise_negative_variance_guard():
  # When teps is small enough that v_e = teps - (dobs*sigmaa/sqrt(scatfrac))^2
  # is <= 0, the target must contribute exactly zero for that reflection
  # rather than producing a NaN/inf or spurious gradient.
  n_refl = 3
  inputs = make_inputs(n_refl, centric=False, seed=5)
  # Force d^2 to exceed teps for every reflection.
  inputs["dobs"] = flex.double([0.95, 0.9, 0.97])
  inputs["sigmaa"] = flex.double([0.9, 0.9, 0.9])
  inputs["scatfrac"] = flex.double([1.0, 1.0, 1.0])
  inputs["teps"] = flex.double([0.5, 0.5, 0.5])
  result = target_work(inputs, compute_gradients=True)
  for t in result.target_per_reflection():
    assert t == 0.0
  # No work reflections should contribute a nonzero gradient either.
  for g in result.gradients_work():
    assert g == 0j

def exercise_scatfrac_sensitivity():
  # Construct f_calc tuned so that d*fc == feff exactly (d = dobs*sigmaa/
  # sqrt(scatfrac)), i.e. the best possible fit for this scatfrac. Moving
  # scatfrac away from that tuned value should make the fit worse
  # (target_work, a minimize-me quantity, increases), since fc no longer
  # matches the (now different) effective d.
  n_refl = 6
  inputs = make_inputs(n_refl, centric=False, seed=2)
  tuned = dict(inputs)
  fc_mag = flex.double([
    inputs["f_eff"][i]
    / (inputs["dobs"][i] * inputs["sigmaa"][i]
       / math.sqrt(inputs["scatfrac"][i]))
    for i in range(n_refl)])
  tuned["f_calc"] = flex.complex_double([
    complex(fc_mag[i] * math.cos(0.2 * i), fc_mag[i] * math.sin(0.2 * i))
    for i in range(n_refl)])
  base = target_work(tuned).target_work()
  mismatched = dict(tuned)
  mismatched["scatfrac"] = inputs["scatfrac"] * 0.3
  mismatched_val = target_work(mismatched).target_work()
  assert mismatched_val > base, (base, mismatched_val)

def exercise_gradient_descent_direction():
  # NOTE: LLGI is a *gain* relative to a Wilson/null prior, so "a
  # particular hand-picked |Fcalc| beats some other hand-picked |Fcalc|"
  # is not a safe premise to test directly -- the maximum-likelihood
  # |Fcalc| for a Rice/Sim-weighted target is not the naive feff/d point
  # estimate (that ignores the Bessel-function bias correction the
  # analytic gradient itself encodes; see mlf.h and llgi.h's use of
  # i1_over_i0/tanh). The property that *is* fundamental and safe to
  # assert directly: stepping f_calc a small amount in the direction the
  # analytic gradient (from gradients_work()) points should *decrease*
  # target_work (since gradients_work() is, by this file's convention,
  # d(target)/d(conj(f_calc)) with target already the minimize-me
  # quantity) -- i.e. the analytic gradient is consistent with target_work
  # for use in gradient-descent-based refinement (this is effectively a
  # coarse, aggregate cross-check of the finite-difference gradient test
  # above, phrased as a descent-direction sanity check rather than a
  # per-component comparison).
  n_refl = 6
  inputs = make_inputs(n_refl, centric=False, seed=4)
  result = target_work(inputs, compute_gradients=True)
  base_val = result.target_work()
  grads = result.gradients_work()
  work_indices = [i for i in range(n_refl) if not inputs["r_free_flags"][i]]
  step = 1.e-4
  stepped = dict(inputs)
  f_calc_stepped = inputs["f_calc"].deep_copy()
  for k, ih in enumerate(work_indices):
    g = grads[k]
    # The finite-difference test above establishes g.real == d(target)/
    # d(Re fc) and g.imag == d(target)/d(Im fc) directly (no extra
    # conjugation needed), so the descent step is simply -step * g.
    f_calc_stepped[ih] = f_calc_stepped[ih] - step * g
  stepped["f_calc"] = f_calc_stepped
  stepped_val = target_work(stepped).target_work()
  assert stepped_val < base_val, (base_val, stepped_val)

def exercise_sigmaa_scatfrac_finite_difference_gradients():
  # Verifies llgi.h's d_target_one_h_over_sigmaa_scatfrac (used by the
  # sigmaA(resolution)/ScatFrac(resolution) B-spline estimator, see
  # doc/llgi_target_design.md sec. 5.2) against central finite differences
  # of llgi_sigmaa_scatfrac_target_and_gradients' summed target -- through
  # the actual compiled/bound code path, not just the standalone hand
  # derivation this function's math was originally checked against.
  # Mixes centric and acentric reflections in one call (unlike
  # exercise_finite_difference_gradients, which uses a uniform flag per
  # call) since d_target_one_h_over_sigmaa_scatfrac's centric/acentric
  # branches share no code with each other or with
  # d_target_one_h_over_fc's branches.
  n_refl = 8
  inputs = make_inputs(n_refl, centric=False, seed=7)
  inputs["centric_flags"] = flex.bool(
    [(i % 3 == 0) for i in range(n_refl)])
  selection = flex.bool([True] * n_refl)

  def run(sigmaa, scatfrac):
    return ext.llgi_sigmaa_scatfrac_target_and_gradients(
      f_eff=inputs["f_eff"],
      selection=selection,
      f_calc=inputs["f_calc"],
      dobs=inputs["dobs"],
      sigmaa=sigmaa,
      scatfrac=scatfrac,
      scale_factor=1.0,
      teps=inputs["teps"],
      resn=inputs["resn"],
      centric_flags=inputs["centric_flags"])

  result = run(inputs["sigmaa"], inputs["scatfrac"])
  ana_dsigmaa = result.d_target_by_dsigmaa()
  ana_dscatfrac = result.d_target_by_dscatfrac()
  eps = 1.e-6
  for i in range(n_refl):
    sa_p = inputs["sigmaa"].deep_copy(); sa_p[i] += eps
    sa_m = inputs["sigmaa"].deep_copy(); sa_m[i] -= eps
    fd_dsigmaa = (run(sa_p, inputs["scatfrac"]).target()
                  - run(sa_m, inputs["scatfrac"]).target()) / (2 * eps)
    assert approx_equal(ana_dsigmaa[i], fd_dsigmaa, eps=5.e-5), (
      i, ana_dsigmaa[i], fd_dsigmaa)
    sf_p = inputs["scatfrac"].deep_copy(); sf_p[i] += eps
    sf_m = inputs["scatfrac"].deep_copy(); sf_m[i] -= eps
    fd_dscatfrac = (run(inputs["sigmaa"], sf_p).target()
                    - run(inputs["sigmaa"], sf_m).target()) / (2 * eps)
    assert approx_equal(ana_dscatfrac[i], fd_dscatfrac, eps=5.e-5), (
      i, ana_dscatfrac[i], fd_dscatfrac)

def _reference_target_original_three_term_form(
      f_eff, f_calc, dobs, sigmaa, scatfrac, k, teps, resn, centric):
  """ The ORIGINAL, phasertng-shaped LLGI target for one reflection: the
  Wilson/null-hypothesis baseline carried as a separate additive term
  wll = EOBS^2/teps = feff^2/(teps^2*resn^2), added after subtracting
  feff^2/V. llgi.h now instead folds those two feff terms together into a
  single -(d^2/teps)*feff^2/V (see its target_one_h comments); this is an
  independent reimplementation of the pre-simplification form, kept here
  so the equivalence is pinned by a test rather than only by the algebra.
  """
  from scitbx.math import bessel_ln_of_i0
  d = dobs * (sigmaa / math.sqrt(scatfrac)) * k
  resn_sq = resn * resn
  v_e = teps - d * d
  if(v_e <= 0.0):
    return 0.0
  v = teps * resn_sq * v_e
  teps_resn_sq = teps * resn_sq
  feff_sq = f_eff * f_eff
  ec_sq = (d * f_calc) ** 2
  x = 2.0 * f_eff * (d * f_calc) / v
  wll = feff_sq / (teps * teps_resn_sq)
  if(not centric):
    ll = -(math.log(v / (teps * teps_resn_sq)) + (feff_sq + ec_sq) / v)
    ll += bessel_ln_of_i0(x)
    ll += wll
  else:
    wll /= 2.0
    ll_core = -(math.log(v / (teps * teps_resn_sq)) + (feff_sq + ec_sq) / v)
    x_half = x / 2.0
    ln_cosh = x_half + math.log((1.0 + math.exp(-2.0 * x_half)) / 2.0)
    ll = ll_core / 2.0 + ln_cosh + wll
  return -ll

def exercise_symmetrised_form_matches_original():
  # llgi.h's target_one_h was simplified to be symmetric in feff and fc:
  # the separate additive wll baseline was folded into the feff^2 term,
  # giving -(d^2/teps)*feff^2/V as the partner of -(d*fc)^2/V. Verify that
  # simplification changes nothing, against an independent Python
  # reimplementation of the original three-term form.
  #
  # The teps != 1 cases matter specifically: on the F-scale the folded
  # multiplier is d^2/TEPS, not the plain d^2 that the E-scale form
  # (llgi_e.h, where TEPS is identically 1) uses. The two agree only at
  # TEPS == 1, so a plain-d^2 F-scale implementation would pass every
  # existing test (all of which use teps == 1) and be wrong precisely
  # where deferred tNCS support is meant to land. At teps = 2 the two
  # differ by more than a factor of two.
  import random
  rnd = random.Random(7)
  n_checked = 0
  for trial in range(400):
    f_eff = rnd.uniform(0.05, 30.0)
    fc_mag = rnd.uniform(0.05, 30.0)
    dobs = rnd.uniform(0.4, 0.99)
    sigmaa = rnd.uniform(0.05, 0.95)
    scatfrac = rnd.uniform(0.3, 1.4)
    k = rnd.uniform(0.5, 1.5)
    # Deliberately include teps != 1 (the case no other test covers).
    teps = 1.0 if (trial % 2 == 0) else rnd.uniform(1.0, 3.0)
    resn = rnd.uniform(0.5, 20.0)
    for centric in (False, True):
      expected = _reference_target_original_three_term_form(
        f_eff, fc_mag, dobs, sigmaa, scatfrac, k, teps, resn, centric)
      phase = 0.37 * trial
      result = ext.llgi_target_and_gradients(
        f_eff         = flex.double([f_eff]),
        r_free_flags  = flex.bool([False]),
        f_calc        = flex.complex_double([complex(
          fc_mag * math.cos(phase), fc_mag * math.sin(phase))]),
        dobs          = flex.double([dobs]),
        sigmaa        = flex.double([sigmaa]),
        scatfrac      = flex.double([scatfrac]),
        scale_factor  = k,
        teps          = flex.double([teps]),
        resn          = flex.double([resn]),
        centric_flags = flex.bool([centric]),
        compute_gradients = False)
      # target() is the sum over the work set (one reflection here).
      assert approx_equal(result.target(), expected, eps=1.e-9), (
        result.target(), expected, teps, centric)
      n_checked += 1
  assert n_checked == 800, n_checked  # 400 trials x {acentric, centric}

def exercise_symmetrised_form_is_symmetric_in_feff_and_fc():
  # The point of the simplification: with the baseline folded in, the
  # target depends on feff and fc through the symmetric pair
  # (d^2/teps)*feff^2 and (d*fc)^2 (plus the Bessel argument, itself
  # symmetric in the two). At teps == 1 and k == 1 the substitution
  # feff <-> fc must therefore leave the target unchanged.
  import random
  rnd = random.Random(13)
  for trial in range(200):
    a = rnd.uniform(0.05, 20.0)
    b = rnd.uniform(0.05, 20.0)
    dobs = rnd.uniform(0.4, 0.99)
    sigmaa = rnd.uniform(0.05, 0.9)
    for centric in (False, True):
      def t(f_eff, fc_mag):
        return ext.llgi_target_and_gradients(
          f_eff         = flex.double([f_eff]),
          r_free_flags  = flex.bool([False]),
          f_calc        = flex.complex_double([complex(fc_mag, 0.0)]),
          dobs          = flex.double([dobs]),
          sigmaa        = flex.double([sigmaa]),
          scatfrac      = flex.double([1.0]),
          scale_factor  = 1.0,
          teps          = flex.double([1.0]),
          resn          = flex.double([1.0]),
          centric_flags = flex.bool([centric]),
          compute_gradients = False).target()
      assert approx_equal(t(a, b), t(b, a), eps=1.e-10), (
        t(a, b), t(b, a), a, b, centric)

def exercise():
  exercise_finite_difference_gradients(centric=False)
  exercise_finite_difference_gradients(centric=True)
  exercise_resn_scale_invariance()
  exercise_negative_variance_guard()
  exercise_scatfrac_sensitivity()
  exercise_gradient_descent_direction()
  exercise_sigmaa_scatfrac_finite_difference_gradients()
  exercise_symmetrised_form_matches_original()
  exercise_symmetrised_form_is_symmetric_in_feff_and_fc()
  print("OK")

if (__name__ == "__main__"):
  exercise()
