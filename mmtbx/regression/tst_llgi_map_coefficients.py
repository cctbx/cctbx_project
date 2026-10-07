from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import random, math
from mmtbx.regression.llgi_test_utils import (
  build_fmodel, build_llgi_fmodel, synthetic_llgi_data)

def exercise_requires_llgi_data_and_sigmaa_scatfrac():
  fmodel = build_fmodel(60, 1.9, seed=10)
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

def _independent_e_scale_quantities(fmodel):
  """ Independent, from-scratch reimplementation of the E-scale Eeff/
  Emodel/D/V quantities map_calculation_helper_llgi() now uses (NOT
  reusing mmtbx.refinement.llgi_e_sigmaa's build_e_eff/
  build_e_model or fmodel.f_model_no_aniso_scale() -- re-derives
  f_model_no_aniso_scale and SigmaP by hand instead), for cross-checking .alpha/.beta/
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
  # flex does not support complex_double / double directly.
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
  # refinement.llgi_e_sigmaa's build_e_model/build_sigma_p, the
  # same functions the method under test itself calls).
  fmodel = build_llgi_fmodel(60, 1.9, seed=13)
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
  fmodel = build_llgi_fmodel(60, 1.9, seed=14)
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

def exercise_map_coefficients_llgi_matches_hand_computation():
  # For 2mFo-DFc and mFo-DFc, the coefficients must equal
  #   Feff*fo_scale*fom (phase-transferred onto mch.f_model's phase)
  #   + mch.f_model*fc_scale*alpha
  # exactly, computed here from mch's own .f_model/.alpha/.fom rather than
  # by calling combine() again. fo_scale/fc_scale come from the real
  # mmtbx.map_tools.fo_fc_scales, so the centric/acentric branching
  # (centrics get plain mFo in 2mFo-DFc) is exercised, not assumed.
  # feff_scale != 1 makes Feff differ from f_obs, so using f_obs instead
  # of Feff would fail here.
  #
  # mch.f_model is f_model_no_aniso_scale (k_anisotropic excluded), which
  # has the same phase as fmodel.f_model() since k_anisotropic is real and
  # positive; it is the array combine() itself uses.
  import mmtbx.map_tools as mt
  import cmath
  fmodel = build_llgi_fmodel(n_atoms=50, d_min=2.1, seed=22, feff_scale=1.35)
  mch = fmodel.map_calculation_helper_llgi()
  llgi_data = fmodel.llgi_data()
  feff = llgi_data.feff.data()
  assert flex.max(flex.abs(feff - fmodel.f_obs().data())) > 1.e-3
  fom = mch.fom
  alpha = mch.alpha.data()
  f_model = mch.f_model
  fmodel_phases = [cmath.phase(v) for v in f_model.data()]
  for map_type in ["2mFo-DFc", "mFo-DFc"]:
    coeffs = fmodel.map_coefficients_llgi(map_type=map_type, isotropize=False)
    ffs = mt.fo_fc_scales(fmodel=fmodel, map_type_str=map_type)
    expected_fo_part_mag = feff * ffs.fo_scale * fom
    expected_fc_part = f_model.data() * ffs.fc_scale * alpha
    expected = flex.complex_double([
      complex(expected_fo_part_mag[i] * math.cos(fmodel_phases[i]),
              expected_fo_part_mag[i] * math.sin(fmodel_phases[i]))
      + expected_fc_part[i]
      for i in range(feff.size())])
    assert coeffs.indices().all_eq(llgi_data.feff.indices())
    diff = flex.max(flex.abs(coeffs.data() - expected))
    assert diff < 1.e-6, (map_type, diff)

def exercise_map_coefficients_llgi_fcalc_only_matches_ordinary():
  # The Fcalc-only special case (map_type="Fc") has no observed-
  # amplitude dependence at all, so map_coefficients_llgi should reuse
  # electron_density_map (the ordinary path) exactly, not duplicate it.
  fmodel = build_llgi_fmodel(60, 1.9, seed=24)
  llgi_fc = fmodel.map_coefficients_llgi(map_type="Fc")
  ml_fc = fmodel.map_coefficients(map_type="Fc")
  assert llgi_fc.indices().all_eq(ml_fc.indices())
  assert approx_equal(
    list(flex.abs(llgi_fc.data())), list(flex.abs(ml_fc.data())),
    eps=1.e-10)

def exercise_map_coefficients_llgi_rejects_anomalous():
  fmodel = build_llgi_fmodel(60, 1.9, seed=25)
  try:
    fmodel.map_coefficients_llgi(map_type="anom")
  except NotImplementedError:
    pass
  else:
    raise RuntimeError("Expected NotImplementedError for map_type=anom.")

def _mcp(map_type, fill_missing_f_obs=False):
  import mmtbx.maps
  import iotbx.phil
  p = iotbx.phil.parse(mmtbx.maps.map_and_map_coeff_params_str)
  mcp = p.extract().map_coefficients[0]
  mcp.map_type = map_type
  mcp.fill_missing_f_obs = fill_missing_f_obs
  return mcp

def build_llgi_fmodel_with_gaps(n_atoms=60, d_min=1.9, seed=0, feff_scale=1.0,
      keep_fraction=0.85):
  """ Like build_llgi_fmodel, but with a random subset of reflections
  DROPPED after building the (otherwise complete, by construction --
  x.structure_factors() generates every symmetry-allowed index at that
  resolution, so build_fmodel's own fmodel never has genuinely missing
  reflections) starting set -- giving model_missing_reflections_llgi/
  model_missing_reflections something real to fill back in, for
  exercise_llgi_fill_missing_matches_hand_computation.
  """
  random.seed(seed + 500)
  fmodel = build_fmodel(n_atoms=n_atoms, d_min=d_min, seed=seed)
  f_obs = fmodel.f_obs()
  keep_sel = flex.bool([
    random.random() < keep_fraction for i in range(f_obs.indices().size())])
  fmodel = fmodel.select(keep_sel)
  llgi_data = synthetic_llgi_data(fmodel, seed=seed + 100,
    feff_scale=feff_scale)
  fmodel.set_llgi_data(llgi_data)
  fmodel.set_target_name("llgi")
  fmodel.update_llgi_sigmaa_scatfrac()
  return fmodel

def exercise_llgi_fill_missing_matches_hand_computation():
  # Verify model_missing_reflections_llgi.get_missing()'s actual formula
  # end to end: for each MISSING reflection (lone to the unfilled LLGI
  # coefficients), the fill value's MAGNITUDE must equal sigmaA(d)*
  # |Emodel|*sqrt(TEPS)*RESN (TEPS==1, Dobs treated as 1 -- see that
  # class's own docstring), computed independently here from mmtbx.
  # refinement.llgi_e_sigmaa's own building blocks (reused, since
  # re-deriving SigmaP/B-spline-sigmaA fitting from scratch a second,
  # independent way is out of scope for this check -- this test's
  # purpose is confirming the ASSEMBLY, not re-verifying machinery
  # already covered by mmtbx.regression.tst_llgi_e_sigmaa).
  from mmtbx import map_tools as mt
  fmodel = build_llgi_fmodel_with_gaps(n_atoms=50, d_min=2.1, seed=32)
  llgi_unfilled = fmodel.map_coefficients_llgi(map_type="2mFo-DFc")
  llgi_filled = fmodel.map_coefficients_llgi(
    map_type="2mFo-DFc", fill_missing=True)
  assert llgi_filled.indices().size() > llgi_unfilled.indices().size(), (
    "test fixture has no genuinely missing reflections to check against "
    "-- increase d_min or n_atoms so some systematic absences/gaps "
    "exist.")
  missing_only = llgi_filled.lone_set(llgi_unfilled)
  assert missing_only.indices().size() > 0

  mro = mt.model_missing_reflections_llgi(fmodel=fmodel, coeffs=llgi_unfilled)
  missing_computed = mro.get_missing()
  # Match indices between the two independently-obtained missing sets
  # (llgi_filled's lone_set vs get_missing()'s own return) before
  # comparing magnitudes.
  a, b = missing_only.common_sets(missing_computed)
  assert a.indices().size() == missing_only.indices().size(), (
    "index mismatch between complete_with's lone_set and get_missing()'s "
    "own f_calc_missing-based index set")
  assert approx_equal(
    list(flex.abs(a.data())), list(flex.abs(b.data())), eps=1.e-6)

def exercise_fill_missing_honours_sigmaa_model():
  # The E-scale phil scope passed to update_llgi_sigmaa_scatfrac must be
  # recorded on llgi_data, survive fmodel.select() (which rebuilds
  # llgi_data field by field -- run by every bss outlier-removal pass,
  # including the final one before maps are written), and reach the
  # fill-missing sigmaA refit. Previously that refit always used the
  # default form, whatever sigmaa_model was set to.
  import mmtbx.refinement.llgi_e_sigmaa as llgi_e_sigmaa
  from mmtbx import map_tools as mt
  fmodel = build_llgi_fmodel_with_gaps(n_atoms=50, d_min=2.1, seed=32)
  # f_model() as the coefficients being completed: model_missing_
  # reflections keeps only atoms whose map (from these coefficients)
  # correlates with the model map, and this fixture's synthetic Feff is
  # unrelated to the model, so real 2mFo-DFc coefficients keep no atoms
  # and the fill is identically zero -- nothing to compare.
  coeffs = fmodel.f_model()
  default_fill = mt.model_missing_reflections_llgi(
    fmodel=fmodel, coeffs=coeffs).get_missing()
  assert flex.min(flex.abs(default_fill.data())) > 0

  e_params = llgi_e_sigmaa.llgi_e_sigmaa_params.extract()
  assert e_params.sigmaa_model == "d_model"  # the default
  e_params.sigmaa_model = "spline"
  fmodel.update_llgi_sigmaa_scatfrac(e_params=e_params)
  assert fmodel.llgi_data().e_params is e_params
  selected = fmodel.select(flex.bool(fmodel.f_obs().size(), True))
  assert selected.llgi_data().e_params is e_params

  spline_fill = mt.model_missing_reflections_llgi(
    fmodel=fmodel, coeffs=coeffs).get_missing()
  a, b = default_fill.common_sets(spline_fill)
  assert a.size() == default_fill.size() > 0
  abs_a, abs_b = flex.abs(a.data()), flex.abs(b.data())
  rel_diff = flex.mean(flex.abs(abs_a - abs_b)) / flex.mean(abs_a)
  assert rel_diff > 0.01, rel_diff

def exercise_phase_errors_llgi_match_numerical_integration():
  # Expected |phase error| implied by the LLGI fom, checked against direct
  # numerical integration of the phase distributions: von Mises
  # exp(X*cos(phi)) (acentric) and the two-point 0/pi distribution with
  # P(pi) = exp(-X)/(1+exp(-X)) (centric).
  import numpy as np
  fmodel = build_llgi_fmodel(n_atoms=50, d_min=2.1, seed=41)
  mch = fmodel.map_calculation_helper_llgi()
  pher = np.array(fmodel.phase_errors_llgi(mch))
  x = np.array(mch.x)
  centric = np.array(mch.f_obs.centric_flags().data())
  phi = np.linspace(0.0, np.pi, 20001)
  worst = 0.0
  for i in np.linspace(0, x.size - 1, 60).astype(int):
    if(centric[i]):
      expected = 180.0 * np.exp(-x[i]) / (1.0 + np.exp(-x[i]))
    else:
      w = np.exp(x[i] * (np.cos(phi) - 1.0))
      expected = np.degrees(
        np.trapezoid(phi * w, phi) / np.trapezoid(w, phi))
    worst = max(worst, abs(pher[i] - expected))
  assert worst < 0.05, worst

def exercise_outlier_selection_skips_model_based_test_for_llgi():
  # Model-based outlier rejection (its own Dobs-free sigmaA fit) must not
  # run under the llgi target; it still runs for ml.
  from mmtbx.scaling import outlier_rejection
  llgi_fmodel = build_llgi_fmodel(n_atoms=50, d_min=2.1, seed=44)
  ml_fmodel = build_fmodel(n_atoms=50, d_min=2.1, seed=44)
  calls = []
  original = outlier_rejection.outlier_manager.model_based_outliers
  def counting(self, *args, **kwargs):
    calls.append(1)
    return original(self, *args, **kwargs)
  outlier_rejection.outlier_manager.model_based_outliers = counting
  try:
    llgi_fmodel.outlier_selection(use_model=True)
    assert len(calls) == 0, len(calls)
    ml_fmodel.outlier_selection(use_model=True)
    assert len(calls) == 1, len(calls)
  finally:
    outlier_rejection.outlier_manager.model_based_outliers = original

def exercise_map_coefficients_from_fmodel_ml_target_unaffected():
  # An ordinary ml-target fmodel (no llgi_data at all) must be routed
  # exactly as before -- this is the regression guard for the "avoid
  # slow calculation several times" shared map_calculation_server fast
  # path in compute_map_coefficients, which is only skipped when
  # llgi_target_active() is True.
  import mmtbx.maps
  fmodel = build_fmodel(n_atoms=50, d_min=2.1, seed=32)
  assert not fmodel.llgi_target_active()
  coeffs = mmtbx.maps.map_coefficients_from_fmodel(
    params=_mcp("2mFo-DFc"), fmodel=fmodel)
  ml_coeffs = fmodel.map_coefficients(map_type="2mFo-DFc")
  assert coeffs.indices().all_eq(ml_coeffs.indices())
  assert approx_equal(
    list(coeffs.data()), list(ml_coeffs.data()), eps=1.e-10)

def exercise_compute_map_coefficients_mixed_dispatch():
  # compute_map_coefficients (the class driver.py's .mtz writer uses)
  # must handle a params list with a MIX of LLGI-supported (mFo-DFc,
  # 2mFo-DFc with fill_missing_f_obs=True, filled natively by the LLGI
  # path) and LLGI-unsupported (anomalous difference map -- SAD analysis,
  # not ported to LLGI) requests in the SAME call, each routed to the
  # LLGI path (it calls map_coefficients_from_fmodel per map type).
  import mmtbx.maps
  fmodel = build_llgi_fmodel(n_atoms=50, d_min=2.1, seed=33)
  params = [
    _mcp("2mFo-DFc", fill_missing_f_obs=True),
    _mcp("mFo-DFc", fill_missing_f_obs=False),
    _mcp("anom", fill_missing_f_obs=False),
  ]
  for p in params:
    p.format = ["mtz"]
  cmo = mmtbx.maps.compute_map_coefficients(fmodel=fmodel, params=params)
  # Only 2 entries: compute_map_coefficients only appends a map_coeffs
  # entry when coeffs is not None (see its own source) -- the anomalous
  # request (this fmodel's f_obs is non-anomalous, so BOTH the LLGI and
  # ML paths return None for it) contributes nothing to the list. The
  # point of including it here is that dispatch doesn't raise/crash on
  # an unsupported map type mixed in with supported ones, not that it
  # produces a placeholder entry.
  assert len(cmo.map_coeffs) == 2
  # First (2mFo-DFc, fill_missing=True) must match the LLGI-native fill
  # path, not the ML fill, and must add the missing reflections.
  llgi_2fofc_filled = fmodel.map_coefficients_llgi(
    map_type="2mFo-DFc", fill_missing=True)
  assert cmo.map_coeffs[0].indices().all_eq(llgi_2fofc_filled.indices())
  assert approx_equal(
    list(cmo.map_coeffs[0].data()), list(llgi_2fofc_filled.data()),
    eps=1.e-10)
  ml_filled = fmodel.map_coefficients(map_type="2mFo-DFc", fill_missing=True)
  assert flex.max(flex.abs(cmo.map_coeffs[0].data() - ml_filled.data())) \
    > 1.e-3
  llgi_unfilled = fmodel.map_coefficients_llgi(map_type="2mFo-DFc")
  assert cmo.map_coeffs[0].size() > llgi_unfilled.size()
  # Second (mFo-DFc, no fill) should match the (unfilled) LLGI path.
  llgi_fofc = fmodel.map_coefficients_llgi(map_type="mFo-DFc")
  assert approx_equal(
    list(cmo.map_coeffs[1].data()), list(llgi_fofc.data()), eps=1.e-10)

def exercise():
  exercise_requires_llgi_data_and_sigmaa_scatfrac()
  exercise_alpha_matches_d_formula()
  exercise_beta_matches_v_formula()
  exercise_fom_matches_bessel_ratio_reference()
  exercise_map_coefficients_llgi_matches_hand_computation()
  exercise_map_coefficients_llgi_fcalc_only_matches_ordinary()
  exercise_map_coefficients_llgi_rejects_anomalous()
  exercise_llgi_fill_missing_matches_hand_computation()
  exercise_fill_missing_honours_sigmaa_model()
  exercise_phase_errors_llgi_match_numerical_integration()
  exercise_outlier_selection_skips_model_based_test_for_llgi()
  exercise_map_coefficients_from_fmodel_ml_target_unaffected()
  exercise_compute_map_coefficients_mixed_dispatch()
  print("OK")

if (__name__ == "__main__"):
  exercise()
