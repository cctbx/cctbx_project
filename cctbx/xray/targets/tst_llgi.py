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

def exercise():
  exercise_finite_difference_gradients(centric=False)
  exercise_finite_difference_gradients(centric=True)
  exercise_resn_scale_invariance()
  exercise_negative_variance_guard()
  exercise_scatfrac_sensitivity()
  exercise_gradient_descent_direction()
  print("OK")

if (__name__ == "__main__"):
  exercise()
