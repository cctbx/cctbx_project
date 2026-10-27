from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import mmtbx.refinement.llgi_e_sigmaa as llgi_e_sigmaa
from libtbx.test_utils import approx_equal
import random
from mmtbx.regression.llgi_test_utils import build_fmodel, synthetic_llgi_data

def exercise_f_model_no_aniso_scale():
  # fmodel.f_model_no_aniso_scale() is k_isotropic*(f_calc + k_mask*f_mask
  # + f_part1 + f_part2), and times k_anisotropic gives f_model().
  fmodel = build_fmodel(40, 1.8, seed=1)
  fmnas = fmodel.f_model_no_aniso_scale()
  assert fmnas.indices().all_eq(fmodel.f_obs().indices())
  bulk = fmodel.f_masks()[0].data() * fmodel.k_masks()[0]
  for j in range(1, len(fmodel.f_masks())):
    bulk = bulk + fmodel.f_masks()[j].data() * fmodel.k_masks()[j]
  fmnas_direct = fmodel.k_isotropic() * (fmodel.f_calc().data() + bulk +
    fmodel.f_part1().data() + fmodel.f_part2().data())
  diff = flex.max(flex.abs(fmnas.data() - fmnas_direct))
  assert diff < 1.e-10, diff
  recovered = fmnas.data() * fmodel.k_anisotropic()
  diff = flex.max(flex.abs(recovered - fmodel.f_model().data()))
  assert diff < 1.e-10, diff

def _e_model_unit_resn(fmodel, model_anisotropy="off"):
  # build_e_model for fmodel with RESN = sqrt(EPS) and TEPS = 1
  f_obs = fmodel.f_obs()
  return llgi_e_sigmaa.build_e_model(fmodel.f_model().data(),
    fmodel.k_anisotropic(), f_obs,
    flex.sqrt(f_obs.epsilons().data().as_double()),
    flex.double(f_obs.size(), 1), model_anisotropy=model_anisotropy)

def exercise_e_model_mean_square_one_by_resolution():
  # with RESN = sqrt(EPS), mean |Emodel|^2 EPS is 1 in each of the moments
  # bins (to the linear interpolation between bin means), and Emodel goes
  # as 1/sqrt(EPS)
  fmodel = build_fmodel(n_atoms=60, d_min=1.6, seed=6)
  f_obs = fmodel.f_obs()
  eps = f_obs.epsilons().data().as_double()
  em = _e_model_unit_resn(fmodel)
  import numpy as np
  x = (flex.norm(em.e_model) * eps).as_numpy_array()
  for b in np.array_split(np.argsort(f_obs.d_star_sq().data().as_numpy_array()),
                          6):
    assert 0.9 < x[b].mean() < 1.1, x[b].mean()
  ratio = flex.abs(em.e_model) * flex.sqrt(eps) / flex.abs(
    fmodel.f_model().data())
  assert approx_equal(ratio * flex.sqrt(em.scatfrac), flex.double(
    ratio.size(), 1), eps=1.e-9)
  assert approx_equal(em.e_per_f, flex.abs(em.e_model) / flex.abs(
    fmodel.f_model_no_aniso_scale().data()), eps=1.e-9)

def exercise_e_model_zero_where_undefined():
  # RESN or TEPS <= 0: Emodel = 0, no failure
  fmodel = build_fmodel(n_atoms=30, d_min=2.5, seed=8)
  f_obs = fmodel.f_obs()
  n = f_obs.size()
  resn = flex.double(n, 1)
  teps = flex.double(n, 1)
  resn[0] = 0
  teps[1] = -1
  em = llgi_e_sigmaa.build_e_model(fmodel.f_model().data(),
    fmodel.k_anisotropic(), f_obs, resn, teps)
  assert abs(em.e_model[0]) == 0 and abs(em.e_model[1]) == 0
  assert flex.min(flex.abs(em.e_model).select(flex.size_t(range(2, n)))) > 0

def _synthetic_llgi_inputs(fmodel, seed=10):
  # (dobs, feff, resn) as flex.double on fmodel.f_obs()'s index set.
  d = synthetic_llgi_data(fmodel, seed=seed)
  return d.dobs.data(), d.feff.data(), d.resn.data()

def exercise_bss_k_sol_b_sol_recovers_known_values():
  # A k_mask with known (k_sol, b_sol) and 1% noise is recovered.
  from mmtbx.f_model import ext as f_model_ext
  fmodel = build_fmodel(n_atoms=30, d_min=2.0, seed=14)
  ss = fmodel.f_obs().sin_theta_over_lambda_sq().data()
  true_k_sol, true_b_sol = 0.35, 45.0
  k_mask_true = f_model_ext.k_mask(ss, true_k_sol, true_b_sol)
  rnd = random.Random(13)
  n = ss.size()
  noise = flex.double([rnd.gauss(0.0, 0.01) for i in range(n)])
  k_mask_noisy = k_mask_true * (flex.double(n, 1.0) + noise)
  fmodel.update(k_mask=[k_mask_noisy])

  k_sol, b_sol = llgi_e_sigmaa.bss_k_sol_b_sol(fmodel)
  assert approx_equal(k_sol, true_k_sol, eps=0.02), (k_sol, true_k_sol)
  assert approx_equal(b_sol, true_b_sol, eps=1.0), (b_sol, true_b_sol)

def exercise_bss_k_sol_b_sol_stays_in_bounds():
  # k_mask ~ 0 (a small structure with little solvent): the estimate stays
  # within [0, 0.6]/[0, 150].
  fmodel = build_fmodel(n_atoms=15, d_min=2.0, seed=15)
  k_sol, b_sol = llgi_e_sigmaa.bss_k_sol_b_sol(
    fmodel, k_sol_default=0.35, b_sol_default=46.0)
  assert k_sol is not None and b_sol is not None
  assert 0.0 <= k_sol <= 0.6
  assert 0.0 <= b_sol <= 150.0

def exercise_sigmaa_fit_leaves_k_mask_untouched():
  # For both sigmaA forms, estimate_e_sigmaa_for_fmodel returns one sigmaA
  # in (0, 1) per reflection and leaves fmodel's k_mask unchanged.
  fmodel = build_fmodel(n_atoms=50, d_min=1.75, seed=30)
  dobs, feff, resn = _synthetic_llgi_inputs(fmodel, seed=31)
  k_mask_before = fmodel.k_masks()[0].deep_copy()
  for sigmaa_model in ["spline", "d_model"]:
    params = llgi_e_sigmaa.llgi_e_sigmaa_params.extract()
    params.sigmaa_model = sigmaa_model
    params.sigmaa_max_iterations = 20
    params.d_model_params.max_iterations = 60
    result = llgi_e_sigmaa.estimate_e_sigmaa_for_fmodel(
      fmodel, dobs=dobs, feff=feff, resn=resn, params=params)
    k_mask_after = flex.double(fmodel.k_masks()[0])
    assert approx_equal(list(k_mask_before), list(k_mask_after), eps=0.0)
    assert result.sigmaa.size() == fmodel.f_obs().size()
    for v in result.sigmaa:
      assert 0.0 < v < 1.0, (sigmaa_model, v)
    # evaluate_at reproduces the fitted curve, and is clamped (spline) or
    # defined (d_model) outside the fitted range.
    d_star_sq = fmodel.f_obs().d_star_sq().data()
    assert approx_equal(result.evaluate_at(d_star_sq), result.sigmaa,
      eps=1.e-12)
    d_out = flex.double([0.5 * flex.min(d_star_sq), 2 * flex.max(d_star_sq)])
    for v in result.evaluate_at(d_out):
      assert 0.0 < v < 1.0, (sigmaa_model, v)

def exercise_spline_curvature_penalty_finite_difference():
  import numpy as np
  random_state = np.random.RandomState(23)
  for trial in range(5):
    n = random_state.randint(3, 10)
    coeffs = random_state.normal(size=n)
    weight = random_state.uniform(0.01, 2.0)
    f0, g0 = llgi_e_sigmaa.spline_curvature_penalty_and_gradient(
      coeffs, weight)
    eps = 1.e-6
    g_fd = np.zeros_like(coeffs)
    for i in range(n):
      cp = coeffs.copy(); cp[i] += eps
      cm = coeffs.copy(); cm[i] -= eps
      fp, _ = llgi_e_sigmaa.spline_curvature_penalty_and_gradient(cp, weight)
      fm, _ = llgi_e_sigmaa.spline_curvature_penalty_and_gradient(cm, weight)
      g_fd[i] = (fp - fm) / (2 * eps)
    assert approx_equal(list(g0), list(g_fd), eps=1.e-6)
  # weight = 0 and fewer than 3 coefficients: no restraint
  f, g = llgi_e_sigmaa.spline_curvature_penalty_and_gradient(
    np.array([1.0, -2.0, 3.0, 0.5]), 0.0)
  assert f == 0.0 and list(g) == [0.0] * 4
  f, g = llgi_e_sigmaa.spline_curvature_penalty_and_gradient(
    np.array([1.0, -2.0]), 5.0)
  assert f == 0.0 and list(g) == [0.0] * 2
  # coefficients linear in the index: zero penalty and gradient
  f, g = llgi_e_sigmaa.spline_curvature_penalty_and_gradient(
    np.linspace(-3.0, 4.0, 8), 10.0)
  assert approx_equal(f, 0.0, eps=1.e-10)
  assert approx_equal(list(g), [0.0] * 8, eps=1.e-10)

def exercise_e_sigmaa_curvature_penalty_gradient_finite_difference():
  # Finite-difference check of e_sigmaa_target_evaluator's gradient with
  # the curvature restraint switched on.
  import numpy as np
  fmodel = build_fmodel(n_atoms=40, d_min=1.8, seed=21)
  dobs, feff, resn = _synthetic_llgi_inputs(fmodel, seed=22)
  f_obs = fmodel.f_obs()
  epsilons = f_obs.epsilons().data().as_double()
  d_star_sq = f_obs.d_star_sq().data()
  centric_flags = f_obs.centric_flags().data()

  e_model = flex.abs(llgi_e_sigmaa.build_e_model(fmodel.f_model().data(),
    fmodel.k_anisotropic(), f_obs, resn, flex.double(f_obs.size(), 1)
    ).e_model)
  e_eff = llgi_e_sigmaa.build_e_eff(feff, resn)
  evaluator = llgi_e_sigmaa.e_sigmaa_target_evaluator(
    e_eff=e_eff, test_selection=fmodel.r_free_flags().data(),
    e_model=e_model, dobs=dobs, centric_flags=centric_flags,
    d_star_sq=d_star_sq, n_sigmaa_coeffs=6,
    max_iterations=0,  # one L-BFGS step; the gradient is checked there
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

def exercise_d_model_sigmaa_is_the_default():
  # sigmaa_model defaults to "d_model", and the plain extract() (used
  # when no params are passed) carries d_model_params.
  params = llgi_e_sigmaa.llgi_e_sigmaa_params.extract()
  assert params.sigmaa_model == "d_model"
  assert params.d_model_params.include_constant_term is True

def exercise_e_model_is_the_target_ecalc():
  # the sigmaA fit and the maps use the F-scale target's Ecalc = k F_model/
  # (RESN sqrt(TEPS ScatFrac)), ScatFrac the moments ScatFrac the fit
  # attaches; e_abs_for gives the same values on the observed reflections
  from mmtbx.regression.llgi_test_utils import build_llgi_fmodel
  import mmtbx.refinement.llgi_scatfrac as llgi_scatfrac
  fmodel = build_llgi_fmodel(40, 2.0, seed=4)
  ld = fmodel.llgi_data()
  f_obs = fmodel.f_obs()
  e_calc = flex.abs(fmodel.f_model().data()) * fmodel.scale_ml_wrapper() / (
    ld.resn.data() * flex.sqrt(ld.teps.data() * ld.scatfrac.data()))
  em = llgi_e_sigmaa.build_e_model_as_fitted(fmodel)
  assert approx_equal(flex.abs(em.e_model), e_calc, eps=1.e-9)
  assert approx_equal(em.scatfrac, ld.scatfrac.data(), eps=1.e-12)
  assert approx_equal(em.e_abs_for(f_obs, fmodel.f_model_no_aniso_scale(
    ).data(), ld.resn.data(), ld.teps.data(), fmodel.k_anisotropic()),
    e_calc, eps=1.e-9)
  # the sigmaA fit used e_calc: refitting reproduces it
  r = llgi_e_sigmaa.estimate_e_sigmaa_for_fmodel(fmodel,
    dobs=ld.dobs.data(), feff=ld.feff.data(), resn=ld.resn.data())
  assert approx_equal(r.e_model, e_calc, eps=1.e-9)
  assert approx_equal(r.sigmaa, ld.sigmaa.data(), eps=1.e-9)
  # scatfrac_moments_model.at reproduces the moments on the same set
  m = llgi_scatfrac.scatfrac_moments_model(fmodel.f_model().data(),
    ld.teps.data(), ld.resn.data(), f_obs.d_star_sq().data(),
    miller_set=f_obs, model_anisotropy="ellipsoid")
  assert approx_equal(m.at(f_obs), m.scatfrac, eps=1.e-9)

def run():
  exercise_f_model_no_aniso_scale()
  exercise_e_model_mean_square_one_by_resolution()
  exercise_e_model_zero_where_undefined()
  exercise_e_model_is_the_target_ecalc()
  exercise_bss_k_sol_b_sol_recovers_known_values()
  exercise_bss_k_sol_b_sol_stays_in_bounds()
  exercise_spline_curvature_penalty_finite_difference()
  exercise_e_sigmaa_curvature_penalty_gradient_finite_difference()
  exercise_sigmaa_fit_leaves_k_mask_untouched()
  exercise_d_model_sigmaa_is_the_default()
  print("OK")

if (__name__ == "__main__"):
  run()
