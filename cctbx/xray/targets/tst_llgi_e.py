from __future__ import absolute_import, division, print_function
from cctbx.xray import ext
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import math

def make_inputs(n_refl, centric, seed=0):
  # Deterministic pseudo-random-ish inputs (no RNG dependency) covering a
  # spread of magnitudes, chosen so that dobs*sigmaa stays comfortably
  # below 1 (i.e. v = 1 - d^2 > 0) for every reflection, avoiding the
  # negative-variance guard in the basic checks. Mirrors tst_llgi.py's
  # make_inputs, but on the E-scale (no fc phase, no resn/scatfrac/teps/k).
  e_eff = flex.double([0.3 + 0.17 * ((seed + i) % 7) for i in range(n_refl)])
  dobs = flex.double([0.4 + 0.05 * ((seed + i) % 6) for i in range(n_refl)])
  sigmaa = flex.double([0.5 + 0.03 * ((seed + i) % 5) for i in range(n_refl)])
  e_model = flex.double(
    [0.4 + 0.11 * ((seed + i) % 5) for i in range(n_refl)])
  centric_flags = flex.bool([centric for i in range(n_refl)])
  selection = flex.bool([True for i in range(n_refl)])
  return dict(
    e_eff=e_eff, selection=selection, e_model=e_model, dobs=dobs,
    sigmaa=sigmaa, centric_flags=centric_flags)

def sigmaa_target(inputs):
  return ext.llgi_e_sigmaa_target_and_gradients(
    e_eff=inputs["e_eff"],
    selection=inputs["selection"],
    e_model=inputs["e_model"],
    dobs=inputs["dobs"],
    sigmaa=inputs["sigmaa"],
    centric_flags=inputs["centric_flags"])

def emodel_target(inputs):
  return ext.llgi_e_emodel_target_and_gradients(
    e_eff=inputs["e_eff"],
    selection=inputs["selection"],
    e_model=inputs["e_model"],
    dobs=inputs["dobs"],
    sigmaa=inputs["sigmaa"],
    centric_flags=inputs["centric_flags"])

def exercise_d_target_over_dsigmaa_finite_difference(centric):
  n_refl = 12
  inputs = make_inputs(n_refl, centric=centric, seed=3)
  result = sigmaa_target(inputs)
  ana_grads = result.d_target_by_dsigmaa()
  eps = 1.e-6
  for i in range(n_refl):
    vals = []
    for signed_eps in (eps, -eps):
      sigmaa_pert = inputs["sigmaa"].deep_copy()
      sigmaa_pert[i] = sigmaa_pert[i] + signed_eps
      inputs_pert = dict(inputs)
      inputs_pert["sigmaa"] = sigmaa_pert
      vals.append(sigmaa_target(inputs_pert).target())
    fin_grad = (vals[0] - vals[1]) / (2 * eps)
    # d_target_by_dsigmaa() is per-reflection (see llgi_e.h docstring),
    # already divided by n_selected -- but target() is also the mean over
    # n_selected, and perturbing sigmaa[i] only affects reflection i's own
    # contribution to that mean, so the per-reflection stored gradient
    # equals d(mean target)/d(sigmaa[i]) directly, matching the finite
    # difference of the mean target above.
    assert approx_equal(ana_grads[i], fin_grad, eps=5.e-5), (
      centric, i, ana_grads[i], fin_grad)

def exercise_d_target_over_demodel_finite_difference(centric):
  n_refl = 12
  inputs = make_inputs(n_refl, centric=centric, seed=5)
  result = emodel_target(inputs)
  ana_grads = result.d_target_by_demodel()
  eps = 1.e-6
  for i in range(n_refl):
    vals = []
    for signed_eps in (eps, -eps):
      e_model_pert = inputs["e_model"].deep_copy()
      e_model_pert[i] = e_model_pert[i] + signed_eps
      inputs_pert = dict(inputs)
      inputs_pert["e_model"] = e_model_pert
      vals.append(emodel_target(inputs_pert).target())
    fin_grad = (vals[0] - vals[1]) / (2 * eps)
    assert approx_equal(ana_grads[i], fin_grad, eps=5.e-5), (
      centric, i, ana_grads[i], fin_grad)

def exercise_sigmaa_and_emodel_targets_agree():
  # Both classes compute the identical target_one_h sum over the same
  # selection; only the gradient differs. Cross-check the target values
  # match exactly for the same inputs.
  inputs = make_inputs(10, centric=False, seed=7)
  t1 = sigmaa_target(inputs).target()
  t2 = emodel_target(inputs).target()
  assert approx_equal(t1, t2, eps=1.e-12)

def exercise_negative_variance_guard():
  # dobs*sigmaa >= 1 => v = 1 - d^2 <= 0 => no contribution (target 0 for
  # that reflection, gradient 0), matching llgi.h's guard.
  inputs = make_inputs(4, centric=False, seed=0)
  sigmaa = inputs["sigmaa"].deep_copy()
  sigmaa[0] = 1.0 / inputs["dobs"][0]  # exactly d == 1
  inputs["sigmaa"] = sigmaa
  result = sigmaa_target(inputs)
  assert result.d_target_by_dsigmaa()[0] == 0.0

def exercise_selection_restricts_reflections():
  # Reflections outside `selection` should not affect target() or
  # gradients at all (relied on by the design's Stage-1/Stage-2 split --
  # R-free-only vs working-set-only -- see design note sec. 5-6).
  inputs = make_inputs(6, centric=False, seed=2)
  selection = flex.bool([True, False, True, False, True, False])
  inputs["selection"] = selection
  result = sigmaa_target(inputs)
  grads = result.d_target_by_dsigmaa()
  for i in range(6):
    if not selection[i]:
      assert grads[i] == 0.0
  # Target should equal the mean over the 3 selected reflections only.
  selected_only = dict(inputs)
  for key in ("e_eff", "e_model", "dobs", "sigmaa", "centric_flags"):
    selected_only[key] = inputs[key].select(selection)
  selected_only["selection"] = flex.bool(
    [True] * selected_only["e_eff"].size())
  t_full = sigmaa_target(inputs).target()
  t_selected = sigmaa_target(selected_only).target()
  assert approx_equal(t_full, t_selected, eps=1.e-10)

def exercise_agrees_with_f_scale_functor_in_degenerate_case():
  # With k=1, teps=1, resn=1, scatfrac=1, feff=Eeff, fc=|Emodel| (real,
  # positive), llgi.h's F-scale target_one_h should reduce to exactly the
  # same value as llgi_e.h's E-scale target_one_h -- both are ports of
  # the same phasertng function.cc, and the F-scale change-of-variables
  # (llgi.h's docstring) is an identity in this degenerate case. This
  # cross-checks the E-scale port against the independently-verified
  # F-scale functor, not just against finite differences of itself.
  for centric in (False, True):
    for seed in range(5):
      eeff = 0.3 + 0.11 * seed
      emodel = 0.4 + 0.09 * seed
      dobs = 0.5 + 0.07 * seed
      sigmaa = 0.3 + 0.05 * seed
      r_e = ext.llgi_e_sigmaa_target_and_gradients(
        e_eff=flex.double([eeff]),
        selection=flex.bool([True]),
        e_model=flex.double([emodel]),
        dobs=flex.double([dobs]),
        sigmaa=flex.double([sigmaa]),
        centric_flags=flex.bool([centric]))
      r_f = ext.llgi_sigmaa_scatfrac_target_and_gradients(
        f_eff=flex.double([eeff]),
        selection=flex.bool([True]),
        f_calc=flex.complex_double([complex(emodel, 0.0)]),
        dobs=flex.double([dobs]),
        sigmaa=flex.double([sigmaa]),
        scatfrac=flex.double([1.0]),
        scale_factor=1.0,
        teps=flex.double([1.0]),
        resn=flex.double([1.0]),
        centric_flags=flex.bool([centric]))
      assert approx_equal(r_e.target(), r_f.target(), eps=1.e-12), (
        centric, seed, r_e.target(), r_f.target())

def exercise_zero_selection_is_safe():
  inputs = make_inputs(3, centric=False, seed=0)
  inputs["selection"] = flex.bool([False, False, False])
  result = sigmaa_target(inputs)
  assert result.target() == 0.0
  for g in result.d_target_by_dsigmaa():
    assert g == 0.0

def exercise():
  exercise_d_target_over_dsigmaa_finite_difference(centric=False)
  exercise_d_target_over_dsigmaa_finite_difference(centric=True)
  exercise_d_target_over_demodel_finite_difference(centric=False)
  exercise_d_target_over_demodel_finite_difference(centric=True)
  exercise_sigmaa_and_emodel_targets_agree()
  exercise_agrees_with_f_scale_functor_in_degenerate_case()
  exercise_negative_variance_guard()
  exercise_selection_restricts_reflections()
  exercise_zero_selection_is_safe()
  print("OK")

if (__name__ == "__main__"):
  exercise()
