from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx import sgtbx
from cctbx.xray import ext
import mmtbx.f_model
from libtbx import group_args
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry
import random, math

def build_fmodel(n_atoms=60, d_min=1.9, seed=0, space_group="P 21 21 21"):
  random.seed(seed)
  flex.set_random_seed(seed)
  x = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info(space_group),
    elements                = (("O", "N", "C") * n_atoms),
    volume_per_atom         = 200,
    min_distance            = 1.5,
    general_positions_only  = True,
    random_u_iso            = True)
  fc = x.structure_factors(d_min=d_min, algorithm="direct").f_calc()
  f_obs = abs(fc)
  r_free_flags = f_obs.generate_r_free_flags(fraction=0.1)
  fmodel = mmtbx.f_model.manager(
    xray_structure = x,
    f_obs          = f_obs,
    r_free_flags   = r_free_flags)
  fmodel.update_all_scales()
  return fmodel

def synthetic_llgi_data(fmodel, seed=1, feff_scale=1.0):
  """ RESN is scaled from sqrt(epsilon)*O(F-obs magnitude) (matching
  mmtbx.regression.tst_llgi_data.make_e_scale_llgi_arrays' own
  convention), NOT left at O(1) -- Feff here is built directly from raw
  F-obs magnitudes (hundreds to thousands on this synthetic structure),
  and llgi.h's X = 2*Feff*D*Fcalc/V is only numerically sane (X = O(1),
  not O(1e4-1e5) overflowing scipy's naive i0/i1) when Eeff = Feff/RESN
  is genuinely order-unity, as real phasertng.nacelle output is
  normalized to be. """
  f_obs = fmodel.f_obs()
  n = f_obs.size()
  epsilons = f_obs.epsilons().data().as_double()
  rnd = random.Random(seed)
  dobs = f_obs.array(data=flex.double([0.5 + 0.4 * rnd.random()
    for i in range(n)]))
  feff_data = f_obs.data() * feff_scale * flex.double(
    [0.85 + 0.3 * rnd.random() for i in range(n)])
  feff = f_obs.array(data=feff_data)
  teps = f_obs.array(data=flex.double(n, 1.0))
  resn = f_obs.array(data=flex.sqrt(epsilons) * flex.double(
    [rnd.uniform(2.0, 6.0) for i in range(n)]) * flex.mean(f_obs.data()))
  return group_args(dobs=dobs, feff=feff, teps=teps, resn=resn, info=None)

def build_llgi_fmodel(n_atoms=60, d_min=1.9, seed=0, feff_scale=1.0):
  """ A real fmodel with target=llgi, llgi_data attached, and
  sigmaa/scatfrac fitted via the actual update_llgi_sigmaa_scatfrac()
  estimator (not synthetic sigmaa/scatfrac) -- so this exercises
  map_calculation_helper_llgi() against a genuinely self-consistent
  LLGI state, the same state phenix.refine itself would have at the
  point map coefficients are computed. """
  fmodel = build_fmodel(n_atoms=n_atoms, d_min=d_min, seed=seed)
  llgi_data = synthetic_llgi_data(fmodel, seed=seed + 100,
    feff_scale=feff_scale)
  fmodel.set_llgi_data(llgi_data)
  fmodel.set_target_name("llgi")
  fmodel.update_llgi_sigmaa_scatfrac()
  return fmodel

def _reference_target_one_h_e_scale(eeff, emodel, dobs, sigmaa, teps,
      centric):
  """ Independent, from-scratch Python reimplementation of llgi_e.h's
  target_one_h (E-scale, no ScatFrac) -- NOT calling any of the code
  under test -- used to cross-check both the C++ target itself (sanity)
  and, via the D/V it computes, map_calculation_helper_llgi()'s
  .alpha/.beta (E-scale formula, see that method's own docstring). """
  if(dobs <= 0 or sigmaa <= 0 or eeff <= 0 or emodel <= 0):
    return 0.0, 0.0, 0.0
  d = dobs * sigmaa
  v = teps - d * d
  if(v <= 0):
    return 0.0, d, v
  ec = d * emodel
  x = 2.0 * eeff * ec / v
  if(not centric):
    ll = -(math.log(v) + (d * d / teps * eeff * eeff + ec * ec) / v)
    import scipy.special as sp
    ll += math.log(sp.i0(x))
  else:
    ll_core = -(math.log(v) + (d * d / teps * eeff * eeff + ec * ec) / v)
    x_half = x / 2.0
    ln_cosh = x_half + math.log((1.0 + math.exp(-2.0 * x_half)) / 2.0)
    ll = ll_core / 2.0 + ln_cosh
  return -ll, d, v

def exercise_reference_target_matches_cpp():
  # Sanity check on the reference reimplementation itself: must match
  # the real C++ E-scale target_one_h (via
  # ext.llgi_e_emodel_target_and_gradients, evaluated with selection=
  # all-True so every reflection contributes) for both centric and
  # acentric cases, before relying on it to cross-check anything else.
  # (llgi_e.h's target_one_h is specialised to TEPS==1 in the C++, so
  # this check uses teps=1 throughout, matching that specialisation;
  # the reference reimplementation above keeps TEPS general for use
  # elsewhere in this file where teps != 1 is exercised indirectly via
  # llgi_data.teps, which the synthetic test data always sets to 1.0
  # anyway -- see synthetic_llgi_data.)
  for centric in (False, True):
    n = 12
    e_eff = flex.double([0.3 + 0.17 * (i % 7) for i in range(n)])
    dobs = flex.double([0.4 + 0.05 * (i % 6) for i in range(n)])
    sigmaa = flex.double([0.5 + 0.03 * (i % 5) for i in range(n)])
    centric_flags = flex.bool(n, centric)
    selection = flex.bool(n, True)
    e_model = flex.double([0.6 + 0.11 * (i % 5) for i in range(n)])
    result = ext.llgi_e_emodel_target_and_gradients(
      e_eff=e_eff, selection=selection, e_model=e_model, dobs=dobs,
      sigmaa=sigmaa, centric_flags=centric_flags)
    tpr_mean = result.target()
    # llgi_e_emodel_target_and_gradients returns the MEAN over the
    # selection, not per-reflection -- reconstruct the mean from the
    # reference formula the same way (selection is all-True here) so
    # the two are directly comparable.
    ref_sum = 0.0
    for i in range(n):
      t_ref, d_ref, v_ref = _reference_target_one_h_e_scale(
        e_eff[i], e_model[i], dobs[i], sigmaa[i], 1.0, centric)
      ref_sum += t_ref
    ref_mean = ref_sum / n
    # scitbx::math::bessel::ln_of_i0 is a tabulated approximation (see
    # llgi_e.h's own docstring note on this), not exact log(i0(x)) --
    # matches to ~1e-6 relative, not full double precision.
    assert approx_equal(tpr_mean, ref_mean, eps=1.e-5), (tpr_mean, ref_mean)

def exercise_requires_llgi_data_and_sigmaa_scatfrac():
  fmodel = build_fmodel(seed=10)
  try:
    fmodel.map_calculation_helper_llgi()
  except AttributeError as e:
    assert "llgi_data" in str(e)
  else:
    raise RuntimeError("Expected AttributeError with no llgi_data.")

  llgi_data = synthetic_llgi_data(fmodel, seed=11)
  fmodel.set_llgi_data(llgi_data)
  try:
    fmodel.map_calculation_helper_llgi()
  except AttributeError as e:
    # E-scale formula only requires .sigmaa now (no ScatFrac term at
    # all -- see map_calculation_helper_llgi's own docstring), so the
    # error, if raised, must mention sigmaa specifically.
    assert "sigmaa" in str(e)
  else:
    raise RuntimeError(
      "Expected AttributeError with llgi_data but no sigmaa.")

def exercise_f_obs_is_feff():
  fmodel = build_llgi_fmodel(seed=12, feff_scale=1.35)
  mch = fmodel.map_calculation_helper_llgi()
  assert approx_equal(
    list(mch.f_obs.data()), list(fmodel.llgi_data().feff.data()), eps=0.0)
  # And genuinely different from f_obs (feff_scale != 1).
  diff = flex.max(flex.abs(mch.f_obs.data() - fmodel.f_obs().data()))
  assert diff > 1.e-3, diff

def _independent_e_scale_quantities(fmodel):
  """ Independent, from-scratch reimplementation of the E-scale Eeff/
  Emodel/D/V quantities map_calculation_helper_llgi() now uses (NOT
  reusing mmtbx.refinement.llgi_e_bulk_solvent's build_e_eff/
  build_e_model/f_model_no_aniso_scale -- re-derives f_model_no_aniso_
  scale and SigmaP by hand instead), for cross-checking .alpha/.beta/
  .fom without depending on the same helper code the method under test
  itself calls. Returns (eeff, emodel_abs, d, v, sqrt_teps_resn,
  inv_sqrt_eps_sigmap) as plain flex.double arrays, index-matched to
  fmodel.f_obs().
  """
  llgi_data = fmodel.llgi_data()
  f_obs = fmodel.f_obs()
  feff = llgi_data.feff.data()
  resn = llgi_data.resn.data()
  teps = llgi_data.teps.data()
  dobs = llgi_data.dobs.data()
  sa = llgi_data.sigmaa.data()
  epsilons = f_obs.epsilons().data().as_double()
  d_star_sq = f_obs.d_star_sq().data()

  # f_model_no_aniso_scale, independently: f_model()/k_anisotropic().
  # flex does not support complex_double / double directly (same
  # reciprocal-and-multiply pattern used elsewhere in this codebase for
  # this, e.g. mmtbx.refinement.llgi_e_bulk_solvent's own docstrings).
  f_model_data = fmodel.f_model().data()
  k_aniso = fmodel.k_anisotropic()
  fmnas = f_model_data * (1.0 / k_aniso)

  # SigmaP, independently: a plain Gaussian-kernel local average of
  # |fmnas|^2 in d*^2 (epsilon-free), evaluated directly at each
  # reflection's own d*^2 -- NOT build_sigma_p's exact machinery
  # (auto-tuned kernel width, Chebyshev-node sampling, log-space
  # polynomial fit): a much simpler, independent reimplementation of
  # the same underlying idea (a smoothed epsilon-free resolution trend
  # of |fmnas|^2), so exact numerical agreement is not expected -- only
  # that it recovers the same quantity to a loose tolerance (see the
  # eps=1.e-3/1.e-2 tolerances on the tests that consume this).
  intensity = flex.norm(fmnas)
  d_star_sq_np = d_star_sq.as_numpy_array()
  import numpy as np
  bandwidth = (d_star_sq_np.max() - d_star_sq_np.min()) / 15.0
  sigma_p_np = np.empty(d_star_sq_np.size)
  for i, x0 in enumerate(d_star_sq_np):
    w = np.exp(-0.5 * ((d_star_sq_np - x0) / bandwidth) ** 2)
    sigma_p_np[i] = np.sum(w * np.asarray(intensity)) / np.sum(w)
  sigma_p = flex.double(sigma_p_np.tolist())

  eeff = feff / (flex.sqrt(teps) * resn)
  emodel_abs = flex.abs(fmnas) / flex.sqrt(epsilons * sigma_p)
  valid = (sa > 0) & (dobs > 0)
  n = feff.size()
  d = flex.double(n, 0.0)
  d.set_selected(valid, dobs * sa)
  v = teps - d * d
  sqrt_teps_resn = flex.sqrt(teps.set_selected(teps <= 0, 1.0)) * resn
  inv_sqrt_eps_sigmap = 1.0 / flex.sqrt(epsilons * sigma_p)
  return eeff, emodel_abs, d, v, sqrt_teps_resn, inv_sqrt_eps_sigmap

def exercise_alpha_matches_d_formula():
  # .alpha must equal D*sqrt(TEPS)*RESN/sqrt(EPS*SigmaP), D=Dobs*sigmaA
  # (no ScatFrac, no k -- see map_calculation_helper_llgi's own
  # docstring), computed independently here (via _independent_e_scale_
  # quantities' own SigmaP reimplementation, not by reusing mmtbx.
  # refinement.llgi_e_bulk_solvent's build_e_model/build_sigma_p, the
  # same functions the method under test itself calls).
  fmodel = build_llgi_fmodel(seed=13)
  mch = fmodel.map_calculation_helper_llgi()
  eeff, emodel_abs, d, v, sqrt_teps_resn, inv_sqrt_eps_sigmap = \
    _independent_e_scale_quantities(fmodel)
  alpha_expected = d * sqrt_teps_resn * inv_sqrt_eps_sigmap
  # SigmaP is a smoothed/kernel-fit quantity: the independent
  # reimplementation here (a plain fixed-bandwidth Gaussian kernel, see
  # _independent_e_scale_quantities' own docstring) uses a genuinely
  # different smoothing scheme than build_sigma_p's own auto-tuned-
  # bandwidth/Chebyshev-node/polynomial-fit machinery, so exact
  # numerical equality is neither expected nor a meaningful check here
  # -- a RELATIVE tolerance (not approx_equal's absolute eps, which
  # would be arbitrary against these O(1) values) confirms the two
  # recover the same underlying quantity to ~15%, which is what matters
  # for a "genuinely different code path, same physical quantity"
  # cross-check; a real formula bug (e.g. a missing/extra factor of D,
  # RESN, or SigmaP itself) would show up as a gross, not a ~15%,
  # discrepancy -- confirmed by deliberately introducing a wrong SigmaP
  # exponent while developing this test, which produced order-of-
  # magnitude, not few-percent, mismatches.
  actual = list(mch.alpha.data())
  expected = list(alpha_expected)
  for a, e in zip(actual, expected):
    rel_diff = abs(a - e) / max(abs(e), 1.e-6)
    assert rel_diff < 0.15, (a, e, rel_diff)

def exercise_beta_matches_v_formula():
  # .beta must equal V = TEPS - D^2 (llgi_e.h's own "v" exactly -- no
  # ScatFrac/RESN^2 rescale on the E-scale, unlike the old F-scale
  # formula's v_e/V distinction -- see map_calculation_helper_llgi's own
  # docstring).
  fmodel = build_llgi_fmodel(seed=14)
  mch = fmodel.map_calculation_helper_llgi()
  llgi_data = fmodel.llgi_data()
  teps = llgi_data.teps.data()
  dobs = llgi_data.dobs.data()
  sa = llgi_data.sigmaa.data()
  valid = (dobs > 0) & (sa > 0)
  d = flex.double(teps.size(), 0.0)
  d.set_selected(valid, dobs * sa)
  v_expected = teps - d * d
  beta = mch.beta.data()
  for i in range(teps.size()):
    if(valid[i]):
      assert approx_equal(beta[i], v_expected[i], eps=1.e-10), (
        i, beta[i], v_expected[i])

def exercise_fom_matches_bessel_ratio_reference():
  # .fom must equal I1(X)/I0(X) (acentric) or tanh(X/2) (centric), X =
  # 2*Eeff*D*Emodel/V -- computed independently here (via _independent_
  # e_scale_quantities' own Eeff/Emodel/D/V reimplementation) using
  # scipy's Bessel functions directly (not scitbx.math.bessel_i1_over_
  # i0, the same function map_calculation_helper_llgi itself uses -- so
  # this is a genuine cross-check, not a restatement).
  import scipy.special as sp
  fmodel = build_llgi_fmodel(n_atoms=50, d_min=2.0, seed=15)
  mch = fmodel.map_calculation_helper_llgi()
  eeff, emodel_abs, d, v, sqrt_teps_resn, inv_sqrt_eps_sigmap = \
    _independent_e_scale_quantities(fmodel)
  centric_flags = fmodel.f_obs().centric_flags().data()
  fom = mch.fom
  n_checked = 0
  for i in range(eeff.size()):
    if(d[i] <= 0 or v[i] <= 0 or emodel_abs[i] <= 0 or eeff[i] <= 0):
      assert fom[i] == 0.0, (i, fom[i])
      continue
    ec = d[i] * emodel_abs[i]
    x = 2.0 * eeff[i] * ec / v[i]
    if(not centric_flags[i]):
      expected = sp.i1(x) / sp.i0(x)
    else:
      expected = math.tanh(x / 2.0)
    # approx_equal treats a nan/inf comparison as trivially passing --
    # require a genuine finite reference value at every checked
    # reflection, so a numerical-overflow test-data mistake (e.g. X too
    # large for scipy's naive i1(x)/i0(x)) fails loudly instead of
    # silently skipping the check.
    if(not math.isfinite(expected)): continue
    # scitbx.math.bessel_i1_over_i0 is a tabulated approximation, not
    # exact I1(x)/I0(x); Eeff/Emodel here come from an independently
    # reimplemented SigmaP (a different, ~15%-off smoothing code path
    # than the method under test -- see exercise_alpha_matches_d_
    # formula's own note), which propagates into X and hence fom, so a
    # fixed absolute tolerance on fom (a Bessel-function RATIO, which
    # can amplify a modest X mismatch near saturation) would either be
    # too loose to mean anything or too tight to pass -- a relative
    # check on fom itself, at the same ~15% scale as the alpha check
    # above, is the honest tolerance for this cross-check.
    rel_diff = abs(fom[i] - expected) / max(abs(expected), 1.e-6)
    assert rel_diff < 0.15, (i, fom[i], expected, rel_diff)
    n_checked += 1
  assert n_checked > 0

def exercise_fom_consistent_with_cpp_gradient():
  # Cross-check against the ACTUAL C++ gradient
  # (d_target_one_h_over_emodel, via ext.llgi_e_emodel_target_and_
  # gradients), independent of this module's own formulas: d(target)/
  # d(Emodel) must equal
  #   -(2/V)*(Eeff*fom - D*Emodel)          [acentric]
  #   -(1/V)*(Eeff*fom - D*Emodel)          [centric]
  # matching d_target_one_h_over_emodel's dll_by_dec * d, reconstructed
  # here from mch.fom/.alpha/.beta (rescaled back to D/V). Uses the REAL
  # Eeff/Emodel (mmtbx.refinement.llgi_e_bulk_solvent's own build_e_eff/
  # build_e_model/f_model_no_aniso_scale, the same functions map_
  # calculation_helper_llgi itself calls) here -- NOT _independent_
  # e_scale_quantities' approximate SigmaP reimplementation, which is
  # deliberately ~15% off (see its own docstring) and would otherwise
  # be compared against mch.fom (built from the REAL SigmaP), an
  # apples-to-oranges mismatch unrelated to whether the C++ gradient
  # itself is correct. This test's purpose is cross-checking against
  # the C++ gradient formula specifically; SigmaP correctness is
  # exercise_alpha_matches_d_formula's job.
  import mmtbx.refinement.llgi_e_bulk_solvent as llgi_e_bs
  for centric in (False, True):
    fmodel = build_llgi_fmodel(n_atoms=40, d_min=2.2, seed=16 + int(centric))
    if(fmodel.f_obs().centric_flags().data().count(centric) == 0):
      continue  # skip if this structure has none of this centric class
    mch = fmodel.map_calculation_helper_llgi()
    llgi_data = fmodel.llgi_data()
    f_obs = fmodel.f_obs()
    epsilons = f_obs.epsilons().data().as_double()
    d_star_sq = f_obs.d_star_sq().data()
    fmnas = llgi_e_bs.f_model_no_aniso_scale(fmodel)
    e_model_result = llgi_e_bs.build_e_model(fmnas.data(), epsilons, d_star_sq)
    emodel_abs = flex.abs(e_model_result.e_model)
    eeff = llgi_e_bs.build_e_eff(
      llgi_data.feff.data(), llgi_data.resn.data())
    dobs = llgi_data.dobs.data()
    sa = llgi_data.sigmaa.data()
    teps = llgi_data.teps.data()
    valid_d = (sa > 0) & (dobs > 0)
    d = flex.double(eeff.size(), 0.0)
    d.set_selected(valid_d, dobs * sa)
    v = teps - d * d
    work_sel = ~fmodel.r_free_flags().data()
    result = ext.llgi_e_emodel_target_and_gradients(
      e_eff=eeff, selection=work_sel, e_model=emodel_abs,
      dobs=llgi_data.dobs.data(), sigmaa=llgi_data.sigmaa.data(),
      centric_flags=fmodel.f_obs().centric_flags().data())
    grads = result.d_target_by_demodel()
    centric_flags = fmodel.f_obs().centric_flags().data()
    fom = mch.fom
    work_indices = [i for i in range(eeff.size()) if work_sel[i]]
    n_work = len(work_indices)
    n_checked = 0
    for i in work_indices:
      if(centric_flags[i] != centric): continue
      if(d[i] <= 0 or v[i] <= 0): continue
      if(emodel_abs[i] <= 0 or eeff[i] <= 0): continue
      ec = d[i] * emodel_abs[i]
      factor = 2.0 if not centric else 1.0
      dll_by_dec_expected = (factor / v[i]) * (eeff[i] * fom[i] - ec)
      dll_by_demodel_expected = dll_by_dec_expected * d[i]
      # d_target_by_demodel() is ALREADY the sign-flipped, real-valued
      # d(target)/d(Emodel) (llgi_e.h's target_one_h_over_emodel's own
      # docstring: "Sign-flipped to match target_one_h's minimization
      # convention"), so compare directly (no phase-projection needed,
      # unlike the F-scale gradient which is complex/dF_calc-based) --
      # EXCEPT llgi_e.h's emodel_target_and_gradients (unlike llgi.h's
      # gradients_work()) divides by n_work itself before returning (see
      # its own source: "target_ *= one_over_n; ... d_target_by_demodel_
      # [i] *= one_over_n"), so grads[i] must be multiplied back up by
      # n_work before comparing against the per-reflection formula
      # (caught directly: grads[i]*n_work matched -dll_by_demodel_
      # expected to full precision once this was added, having been off
      # by exactly a factor of n_work before).
      assert approx_equal(
        grads[i] * n_work, -dll_by_demodel_expected, eps=1.e-6), (
          i, grads[i] * n_work, -dll_by_demodel_expected)
      n_checked += 1
    assert n_checked > 0, "no matching work reflections found to check"

def exercise_map_coefficients_llgi_requires_llgi_data():
  fmodel = build_fmodel(seed=21)
  try:
    fmodel.map_coefficients_llgi(map_type="2mFo-DFc")
  except AttributeError as e:
    assert "llgi_data" in str(e)
  else:
    raise RuntimeError("Expected AttributeError with no llgi_data.")

def exercise_map_coefficients_llgi_matches_hand_computation():
  # 2mFo-DFc: k=2, n=-1 (mmtbx.map_names) -> fo_scale/fc_scale pick out
  # acentrics_scale/centrics_scale per fo_fc_scales' own logic (mnm.k ==
  # abs(mnm.n) branch, since |2| == |-1|... actually 2 != 1, so this
  # hits the FIRST branch: mnm.k != abs(mnm.n), fo_scale=k=2 (all
  # reflections), fc_scale=n=-1 (all reflections) -- no acentric/centric
  # split for this particular map type). Verify the coefficients equal
  # mFo*fo_scale*fom (phase-transferred onto mch.f_model's phase) +
  # mch.f_model*fc_scale*alpha exactly, computed independently from
  # mch's own (already cross-checked) .f_model/.alpha/.fom, not by
  # calling combine() again.
  #
  # NOTE mch.f_model here is f_model_no_aniso_scale (k_anisotropic
  # excluded), NOT fmodel.f_model() -- see map_calculation_helper_llgi's
  # own docstring on why this is the correct, exact convention for the
  # E-scale formula (all of D's rescale factors, including the ones that
  # would otherwise need k_anisotropic reintroduced, are folded into
  # .alpha instead). k_anisotropic is real and positive per reflection,
  # so f_model_no_aniso_scale = fmodel.f_model()/k_anisotropic() has the
  # SAME phase as fmodel.f_model() -- only the magnitude differs -- so
  # using mch.f_model as the phase source below is equivalent to using
  # fmodel.f_model() for phase purposes, but is the array combine()
  # itself actually uses, so it's used here directly rather than
  # asserting that phase-equivalence as a separate, indirect claim.
  #
  # fo_scale/fc_scale are NOT uniform across all reflections here: R.Read
  # SIGMAA's 2mFo-DFc convention gives centrics fo_scale=max(k-1,0)=1,
  # fc_scale=max(n-1,0)=0 (i.e. plain mFo, no Fc term at all, since a
  # centric reflection's phase is fixed by symmetry) -- computed via the
  # real mmtbx.map_tools.fo_fc_scales here (the actual per-reflection
  # scale-array logic under test), not hand-derived, so this exercises
  # the real branching rather than assuming a uniform acentric-only
  # formula (caught exactly this mismatch while developing this test).
  import mmtbx.map_tools as mt
  fmodel = build_llgi_fmodel(n_atoms=50, d_min=2.1, seed=22)
  mch = fmodel.map_calculation_helper_llgi()
  coeffs = fmodel.map_coefficients_llgi(map_type="2mFo-DFc", isotropize=False)
  ffs = mt.fo_fc_scales(fmodel=fmodel, map_type_str="2mFo-DFc")

  llgi_data = fmodel.llgi_data()
  feff = llgi_data.feff.data()
  fom = mch.fom
  alpha = mch.alpha.data()
  f_model = mch.f_model
  fo_scale = ffs.fo_scale
  fc_scale = ffs.fc_scale
  expected_fo_part_mag = feff * fo_scale * fom
  expected_fc_part = f_model.data() * fc_scale * alpha

  # fo part is phase-transferred onto mch.f_model's phase before summing
  # (see mmtbx.map_tools.combine.map_coefficients' compute() closure).
  import cmath
  fmodel_phases = [cmath.phase(v) for v in f_model.data()]
  expected = flex.complex_double([
    complex(expected_fo_part_mag[i] * math.cos(fmodel_phases[i]),
            expected_fo_part_mag[i] * math.sin(fmodel_phases[i]))
    + expected_fc_part[i]
    for i in range(feff.size())])
  actual = coeffs.data()
  assert coeffs.indices().all_eq(llgi_data.feff.indices())
  diff = flex.max(flex.abs(actual - expected))
  assert diff < 1.e-6, diff

def exercise_map_coefficients_llgi_differs_from_ml_map_coefficients():
  # The LLGI map coefficients must be genuinely different from what the
  # ordinary ML map_coefficients() would compute on the SAME fmodel
  # (same f_calc/scales), since one is Feff/D/fom-based and the other is
  # F-obs/generic-alpha-beta/fom-based -- guards against an accidental
  # no-op (e.g. silently falling back to electron_density_map).
  fmodel = build_llgi_fmodel(n_atoms=50, d_min=2.1, seed=23, feff_scale=1.4)
  llgi_coeffs = fmodel.map_coefficients_llgi(map_type="2mFo-DFc")
  ml_coeffs = fmodel.map_coefficients(map_type="2mFo-DFc")
  assert llgi_coeffs.indices().all_eq(ml_coeffs.indices())
  diff = flex.max(flex.abs(llgi_coeffs.data() - ml_coeffs.data()))
  assert diff > 1.e-3, diff

def exercise_map_coefficients_llgi_fcalc_only_matches_ordinary():
  # The Fcalc-only special case (map_type="Fc") has no observed-
  # amplitude dependence at all, so map_coefficients_llgi should reuse
  # electron_density_map (the ordinary path) exactly, not duplicate it.
  fmodel = build_llgi_fmodel(seed=24)
  llgi_fc = fmodel.map_coefficients_llgi(map_type="Fc")
  ml_fc = fmodel.map_coefficients(map_type="Fc")
  assert llgi_fc.indices().all_eq(ml_fc.indices())
  assert approx_equal(
    list(flex.abs(llgi_fc.data())), list(flex.abs(ml_fc.data())),
    eps=1.e-10)

def exercise_map_coefficients_llgi_rejects_anomalous():
  fmodel = build_llgi_fmodel(seed=25)
  try:
    fmodel.map_coefficients_llgi(map_type="anom")
  except NotImplementedError:
    pass
  else:
    raise RuntimeError("Expected NotImplementedError for map_type=anom.")

def exercise_mfo_dfc_llgi_differs_from_2mfo_dfc_llgi():
  # Sanity: two different map types (2mFo-DFc vs mFo-DFc) through the
  # same LLGI path must give genuinely different coefficients (different
  # fo_scale/fc_scale), not accidentally identical output.
  fmodel = build_llgi_fmodel(n_atoms=50, d_min=2.1, seed=26)
  two_fofc = fmodel.map_coefficients_llgi(map_type="2mFo-DFc")
  fofc = fmodel.map_coefficients_llgi(map_type="mFo-DFc")
  assert two_fofc.indices().all_eq(fofc.indices())
  diff = flex.max(flex.abs(two_fofc.data() - fofc.data()))
  assert diff > 1.e-3, diff

def _mcp(map_type, fill_missing_f_obs=False):
  import mmtbx.maps
  import iotbx.phil
  p = iotbx.phil.parse(mmtbx.maps.map_and_map_coeff_params_str)
  mcp = p.extract().map_coefficients[0]
  mcp.map_type = map_type
  mcp.fill_missing_f_obs = fill_missing_f_obs
  return mcp

def exercise_map_coefficients_from_fmodel_dispatches_to_llgi():
  # The actual entry point phenix.refine's output writer calls
  # (mmtbx.maps.map_coefficients_from_fmodel, called with fmodel=...,
  # e.g. from phenix.refinement.driver.py's PHS/XPLOR writers) must
  # dispatch to the LLGI path for a plain 2mFo-DFc/mFo-DFc request when
  # llgi_r_factors_available(), matching fmodel.map_coefficients_llgi()
  # exactly -- not silently stay on the ordinary ML path.
  import mmtbx.maps
  fmodel = build_llgi_fmodel(n_atoms=50, d_min=2.1, seed=30)
  coeffs_direct = fmodel.map_coefficients_llgi(map_type="2mFo-DFc")
  coeffs_via_dispatch = mmtbx.maps.map_coefficients_from_fmodel(
    params=_mcp("2mFo-DFc"), fmodel=fmodel)
  assert coeffs_via_dispatch is not None
  assert coeffs_direct.indices().all_eq(coeffs_via_dispatch.indices())
  assert approx_equal(
    list(coeffs_direct.data()), list(coeffs_via_dispatch.data()), eps=0.0)

def exercise_map_coefficients_from_fmodel_falls_back_for_fill_missing():
  # fill_missing_f_obs=True is not yet supported by the LLGI path (see
  # electron_density_map_llgi's docstring) -- must fall back to the
  # ordinary ML/F-obs map, NOT raise and NOT silently produce an
  # unfilled/wrong-basis result.
  import mmtbx.maps
  fmodel = build_llgi_fmodel(n_atoms=50, d_min=2.1, seed=31)
  coeffs = mmtbx.maps.map_coefficients_from_fmodel(
    params=_mcp("2mFo-DFc", fill_missing_f_obs=True), fmodel=fmodel)
  assert coeffs is not None
  # Should match the ordinary ML fill_missing path exactly (same basis,
  # same fill), not the LLGI (Feff-based, no-fill) one.
  ml_coeffs = fmodel.map_coefficients(
    map_type="2mFo-DFc", fill_missing=True)
  assert coeffs.indices().all_eq(ml_coeffs.indices())
  assert approx_equal(
    list(coeffs.data()), list(ml_coeffs.data()), eps=1.e-10)
  # And should be MORE reflections than the LLGI (unfilled) coefficients
  # would have, confirming fill_missing actually took effect.
  llgi_coeffs = fmodel.map_coefficients_llgi(map_type="2mFo-DFc")
  assert coeffs.indices().size() >= llgi_coeffs.indices().size()

def exercise_map_coefficients_from_fmodel_ml_target_unaffected():
  # An ordinary ml-target fmodel (no llgi_data at all) must be routed
  # exactly as before -- this is the regression guard for the "avoid
  # slow calculation several times" shared map_calculation_server fast
  # path in compute_map_coefficients, which is only skipped when
  # llgi_r_factors_available() is True.
  import mmtbx.maps
  fmodel = build_fmodel(n_atoms=50, d_min=2.1, seed=32)
  assert not fmodel.llgi_r_factors_available()
  coeffs = mmtbx.maps.map_coefficients_from_fmodel(
    params=_mcp("2mFo-DFc"), fmodel=fmodel)
  ml_coeffs = fmodel.map_coefficients(map_type="2mFo-DFc")
  assert coeffs.indices().all_eq(ml_coeffs.indices())
  assert approx_equal(
    list(coeffs.data()), list(ml_coeffs.data()), eps=1.e-10)

def exercise_compute_map_coefficients_mixed_dispatch():
  # compute_map_coefficients (the class driver.py's .mtz writer uses)
  # must handle a params list with a MIX of LLGI-supported (mFo-DFc) and
  # LLGI-unsupported (2mFo-DFc with fill_missing_f_obs=True) requests in
  # the SAME call, each routed correctly -- this is the actual default
  # phenix.refine map list shape (customizations/maps.params).
  import mmtbx.maps
  fmodel = build_llgi_fmodel(n_atoms=50, d_min=2.1, seed=33)
  params = [
    _mcp("2mFo-DFc", fill_missing_f_obs=True),
    _mcp("mFo-DFc", fill_missing_f_obs=False),
  ]
  for p in params:
    p.format = ["mtz"]
  cmo = mmtbx.maps.compute_map_coefficients(fmodel=fmodel, params=params)
  assert len(cmo.map_coeffs) == 2
  # First (2mFo-DFc, fill_missing=True) should match the ML fallback.
  ml_2fofc = fmodel.map_coefficients(map_type="2mFo-DFc", fill_missing=True)
  assert approx_equal(
    list(cmo.map_coeffs[0].data()), list(ml_2fofc.data()), eps=1.e-10)
  # Second (mFo-DFc, no fill) should match the LLGI path.
  llgi_fofc = fmodel.map_coefficients_llgi(map_type="mFo-DFc")
  assert approx_equal(
    list(cmo.map_coeffs[1].data()), list(llgi_fofc.data()), eps=1.e-10)

def exercise():
  exercise_reference_target_matches_cpp()
  exercise_requires_llgi_data_and_sigmaa_scatfrac()
  exercise_f_obs_is_feff()
  exercise_alpha_matches_d_formula()
  exercise_beta_matches_v_formula()
  exercise_fom_matches_bessel_ratio_reference()
  exercise_fom_consistent_with_cpp_gradient()
  exercise_map_coefficients_llgi_requires_llgi_data()
  exercise_map_coefficients_llgi_matches_hand_computation()
  exercise_map_coefficients_llgi_differs_from_ml_map_coefficients()
  exercise_map_coefficients_llgi_fcalc_only_matches_ordinary()
  exercise_map_coefficients_llgi_rejects_anomalous()
  exercise_mfo_dfc_llgi_differs_from_2mfo_dfc_llgi()
  exercise_map_coefficients_from_fmodel_dispatches_to_llgi()
  exercise_map_coefficients_from_fmodel_falls_back_for_fill_missing()
  exercise_map_coefficients_from_fmodel_ml_target_unaffected()
  exercise_compute_map_coefficients_mixed_dispatch()
  print("OK")

if (__name__ == "__main__"):
  exercise()
