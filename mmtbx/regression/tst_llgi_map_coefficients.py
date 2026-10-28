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
    assert "sigmaa" in str(e)
  else:
    raise RuntimeError(
      "Expected AttributeError with llgi_data but no sigmaa.")

def _independent_e_scale_quantities(fmodel):
  """ The E-scale Eeff, |Emodel|, D, V and scale factors computed without
  the code under test (no build_e_eff/build_e_model or
  fmodel.f_model_no_aniso_scale()), from the target's normalisation (the
  ScatFrac attached by the fit), as flex.double arrays on fmodel.f_obs()'s
  index set: (eeff, emodel_abs, d, v, sqrt_teps_resn, e_per_f), e_per_f =
  |Emodel|/|F_model/k_anisotropic|.
  """
  llgi_data = fmodel.llgi_data()
  feff = llgi_data.feff.data()
  resn = llgi_data.resn.data()
  teps = llgi_data.teps.data()
  dobs = llgi_data.dobs.data()
  sa = llgi_data.sigmaa.data()
  scatfrac = llgi_data.scatfrac.data()
  k = fmodel.scale_ml_wrapper()
  k_aniso = fmodel.k_anisotropic()
  # the target's Ecalc = k F_model/(RESN sqrt(TEPS ScatFrac))
  e_per_f = k * k_aniso / (resn * flex.sqrt(teps * scatfrac))
  eeff = feff / (flex.sqrt(teps) * resn)
  emodel_abs = flex.abs(fmodel.f_model().data()) * k / (
    resn * flex.sqrt(teps * scatfrac))
  valid = (sa > 0) & (dobs > 0)
  n = feff.size()
  d = flex.double(n, 0.0)
  d.set_selected(valid, dobs * sa)
  v = teps - d * d
  sqrt_teps_resn = flex.sqrt(
    teps.deep_copy().set_selected(teps <= 0, 1.0)) * resn
  return eeff, emodel_abs, d, v, sqrt_teps_resn, e_per_f

def exercise_alpha_matches_d_formula():
  # .alpha = D*sqrt(TEPS)*RESN*|Emodel|/|F_model/k_anisotropic|, D =
  # Dobs*sigmaA, against the independent quantities above.
  fmodel = build_llgi_fmodel(60, 1.9, seed=13)
  mch = fmodel.map_calculation_helper_llgi()
  eeff, emodel_abs, d, v, sqrt_teps_resn, e_per_f = \
    _independent_e_scale_quantities(fmodel)
  alpha_expected = d * sqrt_teps_resn * e_per_f
  assert approx_equal(mch.alpha.data(), alpha_expected, eps=1.e-9)

def exercise_beta_matches_v_formula():
  # .beta = V = TEPS - D^2.
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
  # .fom = I1(X)/I0(X) (acentric) or tanh(X/2) (centric), X =
  # 2*Eeff*D*Emodel/V, with scipy's Bessel functions and the independent
  # quantities above.
  import scipy.special as sp
  fmodel = build_llgi_fmodel(n_atoms=50, d_min=2.0, seed=15)
  mch = fmodel.map_calculation_helper_llgi()
  eeff, emodel_abs, d, v, sqrt_teps_resn, e_per_f = \
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
    if(not math.isfinite(expected)): continue
    # scitbx's bessel_i1_over_i0 is an approximation (about 1e-6)
    assert approx_equal(fom[i], expected, eps=1.e-5), (i, fom[i], expected)
    n_checked += 1
  assert n_checked > 0

def exercise_map_coefficients_llgi_matches_hand_computation():
  # 2mFo-DFc and mFo-DFc equal Feff*fo_scale*fom along the model phase
  # plus mch.f_model*fc_scale*alpha (fo_fc_scales gives centrics plain
  # mFo in 2mFo-DFc). feff_scale != 1, so using f_obs instead of Feff
  # would fail. mch.f_model is f_model_no_aniso_scale, with the phase of
  # f_model.
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
  # The Fcalc-only map does not depend on the data: same as the ordinary one.
  fmodel = build_llgi_fmodel(60, 1.9, seed=24)
  llgi_fc = fmodel.map_coefficients_llgi(map_type="Fc")
  ml_fc = fmodel.map_coefficients(map_type="Fc")
  assert llgi_fc.indices().all_eq(ml_fc.indices())
  assert approx_equal(
    list(flex.abs(llgi_fc.data())), list(flex.abs(ml_fc.data())),
    eps=1.e-10)

def exercise_map_coefficients_llgi_anomalous_as_ordinary():
  # Anomalous maps come from F_obs anomalous differences, as without LLGI
  # (here None: the data are not anomalous).
  fmodel = build_llgi_fmodel(60, 1.9, seed=25)
  assert fmodel.map_coefficients_llgi(map_type="anom") is None
  assert fmodel.map_coefficients(map_type="anom") is None

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
  """ build_llgi_fmodel with a random keep_fraction of the reflections, so
  that there are missing reflections to fill. """
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
  # The filled 2mFo-DFc adds exactly model_missing_reflections_llgi's
  # values for the missing reflections.
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
  a, b = missing_only.common_sets(missing_computed)
  assert a.indices().size() == missing_only.indices().size(), (
    "index mismatch between complete_with's lone_set and get_missing()'s "
    "own f_calc_missing-based index set")
  assert approx_equal(
    list(flex.abs(a.data())), list(flex.abs(b.data())), eps=1.e-6)

def exercise_fill_missing_uses_fitted_sigmaa_curve():
  # The fitted sigmaA curve is stored on llgi_data, reproduces llgi_data.
  # sigmaa, survives fmodel.select() (every bss outlier-removal pass) and
  # deep copy/pickling, and gives the fill's sigmaA. The fill follows the
  # sigmaa_model used in the fit.
  import copy, pickle
  import mmtbx.refinement.llgi_e_sigmaa as llgi_e_sigmaa
  from mmtbx import map_tools as mt
  fmodel = build_llgi_fmodel_with_gaps(n_atoms=50, d_min=2.1, seed=32)
  fills = []
  for sigmaa_model in ["d_model", "spline"]:
    e_params = llgi_e_sigmaa.llgi_e_sigmaa_params.extract()
    e_params.sigmaa_model = sigmaa_model
    fmodel.update_llgi_sigmaa_scatfrac(e_params=e_params)
    llgi_data = fmodel.llgi_data()
    curve = llgi_data.sigmaa_curve
    d_star_sq = fmodel.f_obs().d_star_sq().data()
    assert approx_equal(curve(d_star_sq), llgi_data.sigmaa.data(), eps=1.e-12)
    selected = fmodel.select(flex.bool(fmodel.f_obs().size(), True))
    assert selected.llgi_data().sigmaa_curve is curve
    for c in (copy.deepcopy(curve), pickle.loads(pickle.dumps(curve))):
      assert approx_equal(c(d_star_sq), curve(d_star_sq), eps=0)
    mro = mt.model_missing_reflections_llgi(fmodel=fmodel)
    fill = mro.get_missing()
    assert flex.min(flex.abs(fill.data())) > 0
    assert approx_equal(mro.e_scale_missing.sigmaa,
      curve(mro.e_scale_missing.miller_set.d_star_sq().data()), eps=0)
    fills.append(fill)
  a, b = fills[0].common_sets(fills[1])
  assert a.size() == fills[0].size() > 0
  abs_a, abs_b = flex.abs(a.data()), flex.abs(b.data())
  rel_diff = flex.mean(flex.abs(abs_a - abs_b)) / flex.mean(abs_a)
  assert rel_diff > 0.01, rel_diff

def exercise_llgi_fill_is_d_fc():
  # The LLGI fill, evaluated for the observed reflections, is the D*Fc term
  # of the LLGI map coefficients there (mFo - (mFo-DFc), phenix's acentric
  # difference coefficients being doubled), with and without isotropize:
  # the map's own scale, no rescaling. With intensities (the exact
  # posterior: coefficient sigmaA, as for a reflection without data; the
  # Rice form's is Dobs*sigmaA).
  from mmtbx import map_tools as mt
  fmodel = build_llgi_fmodel(n_atoms=40, d_min=2.2, seed=3, rice_kappa=0.1)
  emap = fmodel.electron_density_map_llgi()
  for isotropize in (True, False):
    mfo = emap.map_coefficients(map_type="mFo", isotropize=isotropize)
    diff = emap.map_coefficients(map_type="mFo-DFc", isotropize=isotropize)
    half = flex.double([1.0 if c else 0.5
      for c in mfo.centric_flags().data()])
    dfc = mfo.customized_copy(data=mfo.data() - diff.data() * half)
    fill = mt.model_missing_reflections_llgi(fmodel=fmodel).get_missing(
      isotropize=isotropize)
    f, d = fill.common_sets(dfc)
    assert f.size() == dfc.size()
    assert flex.max(flex.abs(f.data() - d.data()) / flex.abs(d.data())) < 0.01
  # the filled (isotropized) map adds the fill for the reflections it lacks
  fill = mt.model_missing_reflections_llgi(fmodel=fmodel).get_missing()
  filled = emap.map_coefficients(map_type="2mFo-DFc",
    exclude_free_r_reflections=True, fill_missing=True)
  free = fmodel.f_obs().select(fmodel.r_free_flags().data())
  a, b = filled.common_sets(fill.common_set(free))
  assert a.size() == free.size()
  assert approx_equal(a.data(), b.data(), eps=1e-6)
  # times k + n for k*mFo + n*DFc: a difference map is not filled
  work = ~fmodel.r_free_flags().data()
  diff = emap.map_coefficients(map_type="mFo-DFc",
    exclude_free_r_reflections=True, fill_missing=True)
  assert diff.indices().all_eq(fmodel.f_obs().indices().select(work))
  # with TEPS != 1 (no tNCS treatment), nacelle_all on another scale (as
  # read, before apply_scale_k1_to_f_obs) and ANISOBETA: still D*Fc
  import mmtbx.refinement.llgi_hybrid as llgi_hybrid
  from libtbx import group_args
  ld = fmodel.llgi_data()
  f_obs = fmodel.f_obs()
  rnd = random.Random(7)
  teps = f_obs.array(data=flex.double([rnd.uniform(0.5, 2) for i in
    range(f_obs.size())]))
  quad = lambda h: math.exp(-0.002*h[0]**2 + 0.001*h[1]**2 + 0.0005*h[2]**2)
  ab = f_obs.array(data=flex.double([quad(h) for h in f_obs.indices()]))
  fmodel.set_llgi_data(llgi_hybrid.replace_llgi_data(ld, teps=teps,
    anisobeta=ab, nacelle_all=group_args(resn=ld.resn.customized_copy(
      data=ld.resn.data() * 10), anisobeta=ab)))
  emap = fmodel.electron_density_map_llgi()
  mfo = emap.map_coefficients(map_type="mFo")
  diff = emap.map_coefficients(map_type="mFo-DFc")
  half = flex.double([1.0 if c else 0.5 for c in mfo.centric_flags().data()])
  dfc = mfo.customized_copy(data=mfo.data() - diff.data() * half)
  f, d = mt.model_missing_reflections_llgi(fmodel=fmodel).get_missing(
    ).common_sets(dfc)
  assert f.size() == dfc.size()
  assert flex.max(flex.abs(f.data() - d.data()) / flex.abs(d.data())) < 0.01

def exercise_phase_errors_llgi_match_numerical_integration():
  # Rice form (no intensities): expected |phase error| of the fom's phase
  # distribution, checked against direct numerical integration: von Mises
  # exp(X*cos(phi)) (acentric) and the two-point 0/pi distribution with
  # P(pi) = exp(-X)/(1+exp(-X)) (centric).
  import numpy as np
  fmodel = build_llgi_fmodel(n_atoms=50, d_min=2.1, seed=41)
  mch = fmodel.map_calculation_helper_llgi()
  assert mch.phase_error is None and mch.n_exact == 0
  pher = np.array(fmodel.phase_errors_llgi(mch))
  x = np.array(mch.x)
  centric = np.array(mch.f_obs.centric_flags().data())
  for i in np.linspace(0, x.size - 1, 60).astype(int):
    if(centric[i]):
      expected = 180.0 * np.exp(-x[i]) / (1.0 + np.exp(-x[i]))
    else:
      expected = _von_mises_mean_abs_phase(x[i])
    assert abs(pher[i] - expected) < 1.e-6, (i, x[i], pher[i], expected)

def _von_mises_mean_abs_phase(x):
  """ <|phi|> in degrees for p(phi) ~ exp(x cos(phi)) on [-pi, pi], by the
  trapezoid rule on a fine grid over the peak """
  import numpy as np
  upper = min(np.pi, 40. / np.sqrt(x)) if x > 0 else np.pi
  phi = np.linspace(0.0, upper, 200001)
  w = np.exp(x * (np.cos(phi) - 1.0))
  return np.degrees(np.trapezoid(phi * w, phi) / np.trapezoid(w, phi))

def exercise_phase_errors_llgi_exact_posterior_mean():
  # With intensities, the exact posterior: the expected |phase error| is
  # the mean over the posterior of J = |E|^2 of the phase error given J
  # (von Mises with x = kappa sqrt(J), or P(pi) = 1/(1 + exp(2 kappa
  # sqrt(J))) centric), here by the trapezoid rule in s = sqrt(J) and
  # max_lik's von Mises mean (checked in tst_maxlik); larger than the
  # phase error implied by the fom (an |E|-weighted mean cosine).
  import numpy as np
  import scitbx.math
  import mmtbx.refinement.llgi_hybrid as llgi_hybrid
  import mmtbx.refinement.llgi_e_sigmaa as llgi_e_sigmaa
  from mmtbx import max_lik
  fmodel = build_llgi_fmodel(n_atoms=50, d_min=2.1, seed=41, rice_kappa=0.0)
  mch = fmodel.map_calculation_helper_llgi(phase_errors=True)
  assert mch.n_exact > 0
  pher = np.array(fmodel.phase_errors_llgi(mch))
  assert np.all(pher == np.array(mch.phase_error))
  ld = fmodel.llgi_data()
  h = llgi_hybrid.get_exact_data(ld)
  eo, sig = np.array(h.e_obs_sq), np.array(h.sig_e_obs_sq)
  ec = np.array(flex.abs(
    llgi_e_sigmaa.build_e_model_as_fitted(fmodel).e_model))
  sa = np.array(ld.sigmaa.data())
  centric = np.array(mch.f_obs.centric_flags().data())
  fom = np.array(mch.fom)
  exact = np.flatnonzero((sig > 0) & (sa > 0) & (sa < 0.999))
  assert exact.size == mch.n_exact
  n_lower = 0
  for i in exact[np.linspace(0, exact.size - 1, 40).astype(int)]:
    S = 1. - sa[i] ** 2
    b = 0.5 if centric[i] else 1.
    kappa = 2. * b * sa[i] * ec[i] / S
    j_hi = max(eo[i] + 12. * sig[i], 0.)
    j_hi = min(j_hi, (sa[i] * ec[i] + 9. * math.sqrt(S)) ** 2) or j_hi
    s = np.linspace(math.sqrt(max(eo[i] - 12. * sig[i], 0.)),
      math.sqrt(j_hi), 20001)
    j = s * s
    x = kappa * s
    lw = -(j - eo[i]) ** 2 / (2. * sig[i] ** 2) - b * j / S
    if(centric[i]):
      lw += np.logaddexp(x, -x)       # ln cosh, J^(-1/2) dJ = 2 ds
      pe = 180. / (1. + np.exp(np.minimum(2. * x, 700.)))
    else:
      big = x >= 500.
      ln_i0 = np.where(big, x - 0.5 * np.log(2. * np.pi * np.maximum(x, 1.))
        + 1. / (8. * np.maximum(x, 1.)), np.log(np.i0(np.minimum(x, 500.))))
      lw += ln_i0 + np.log(np.maximum(s, 1.e-300))   # dJ = 2 s ds
      one = flex.double(x.size, 1.)
      pe = np.array(max_lik.fom_and_phase_error(f_obs=flex.double(x * 0.5),
        f_model=one, alpha=one, beta=one, epsilons=one,
        centric_flags=flex.bool(x.size, False)).phase_error())
    w = np.exp(lw - lw.max())
    expected = np.trapezoid(w * pe, s) / np.trapezoid(w, s)
    assert abs(pher[i] - expected) < 1.e-5, (i, pher[i], expected)
    # the phase error the fom alone would imply
    if(not centric[i] and fom[i] < 0.999999):
      x_m = scitbx.math.inverse_bessel_i1_over_i0(flex.double([fom[i]]))
      if(pher[i] > _von_mises_mean_abs_phase(x_m[0])): n_lower += 1
  assert n_lower > 0

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
  # An ml-target fmodel gets the ordinary maps.
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
  # compute_map_coefficients (phenix.refine's MTZ output) with the llgi
  # target: filled 2mFo-DFc and mFo-DFc are the LLGI maps, and an
  # anomalous map request (None here, as the data are not anomalous) is
  # handled alongside them.
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
  assert len(cmo.map_coeffs) == 2
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
  exercise_map_coefficients_llgi_anomalous_as_ordinary()
  exercise_llgi_fill_missing_matches_hand_computation()
  exercise_fill_missing_uses_fitted_sigmaa_curve()
  exercise_llgi_fill_is_d_fc()
  exercise_phase_errors_llgi_match_numerical_integration()
  exercise_phase_errors_llgi_exact_posterior_mean()
  exercise_outlier_selection_skips_model_based_test_for_llgi()
  exercise_map_coefficients_from_fmodel_ml_target_unaffected()
  exercise_compute_map_coefficients_mixed_dispatch()
  print("OK")

if (__name__ == "__main__"):
  exercise()
