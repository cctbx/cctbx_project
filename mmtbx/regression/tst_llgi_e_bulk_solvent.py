from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import mmtbx.refinement.llgi_e_bulk_solvent as llgi_e_bs
from libtbx.test_utils import approx_equal
import random
from mmtbx.regression.llgi_test_utils import build_fmodel, synthetic_llgi_data

def exercise_f_model_no_aniso_scale_matches_direct_computation():
  # Reconstruction (f_model()/k_anisotropic()) must exactly match the
  # direct formula k_isotropic*(f_calc + k_mask*f_mask + f_part1 +
  # f_part2), since that is what the C++ core itself computes internally
  # (mmtbx/f_model/f_model.h) before multiplying by k_anisotropic.
  fmodel = build_fmodel(40, 1.8, seed=1)
  fmnas = llgi_e_bs.f_model_no_aniso_scale(fmodel)

  k_iso = fmodel.k_isotropic()
  f_calc = fmodel.f_calc()
  f_masks = fmodel.f_masks()
  k_masks = fmodel.k_masks()
  f_part1 = fmodel.f_part1()
  f_part2 = fmodel.f_part2()
  bulk = f_masks[0].data() * k_masks[0]
  for j in range(1, len(f_masks)):
    bulk = bulk + f_masks[j].data() * k_masks[j]
  fmnas_direct = k_iso * (
    f_calc.data() + bulk + f_part1.data() + f_part2.data())

  diff = flex.max(flex.abs(fmnas.data() - fmnas_direct))
  assert diff < 1.e-10, diff

def exercise_f_model_no_aniso_scale_times_k_aniso_is_f_model():
  # Round-trip sanity check independent of the direct-formula comparison
  # above: multiplying back by k_anisotropic must exactly recover
  # f_model().
  fmodel = build_fmodel(40, 1.8, seed=2)
  fmnas = llgi_e_bs.f_model_no_aniso_scale(fmodel)
  k_aniso = fmodel.k_anisotropic()
  recovered = fmnas.data() * k_aniso
  diff = flex.max(flex.abs(recovered - fmodel.f_model().data()))
  assert diff < 1.e-10, diff

def exercise_sigma_p_constant_intensity_recovers_constant():
  # If |f_model_no_aniso_scale|^2 is (deterministically) constant across
  # all reflections, the smoothed SigmaP trend should recover that
  # constant everywhere -- the simplest possible correctness check on the
  # Chebyshev-node kernel fit, independent of any epsilon handling.
  n = 300
  rnd = random.Random(3)
  d_star_sq = flex.double(sorted(
    rnd.uniform(0.01, 0.3) for i in range(n)))
  const_amplitude = 5.0
  fmnas = flex.complex_double([complex(const_amplitude, 0.0)] * n)
  sigma_p = llgi_e_bs.build_sigma_p(fmnas, d_star_sq, n_nodes=10)
  expected = const_amplitude ** 2
  for v in sigma_p:
    assert approx_equal(v, expected, eps=1.e-6), (v, expected)

def exercise_sigma_p_is_epsilon_free():
  # build_sigma_p must NOT be sensitive to epsilon -- i.e. calling it on
  # the same |fmnas|^2 values regardless of what epsilons the caller will
  # later apply in build_e_model. This is the whole point of calling the
  # low-level kernel_normalisation extension with epsilon forced to
  # all-ones rather than going through the epsilon-dividing Python class
  # (see build_sigma_p's docstring for the earlier draft that got this
  # wrong). Verify indirectly: build_e_model with epsilons=2 everywhere
  # should give exactly the same SigmaP as epsilons=1 everywhere (only
  # the final Emodel division should differ, by a factor of sqrt(2)).
  n = 200
  rnd_x = random.Random(4)
  d_star_sq = flex.double(sorted(
    rnd_x.uniform(0.01, 0.3) for i in range(n)))
  rnd = random.Random(5)
  fmnas = flex.complex_double(
    [complex(rnd.uniform(1.0, 10.0), rnd.uniform(-2.0, 2.0))
     for i in range(n)])
  eps_one = flex.double(n, 1.0)
  eps_two = flex.double(n, 2.0)
  r1 = llgi_e_bs.build_e_model(fmnas, eps_one, d_star_sq)
  r2 = llgi_e_bs.build_e_model(fmnas, eps_two, d_star_sq)
  diff_sigma_p = flex.max(flex.abs(r1.sigma_p - r2.sigma_p))
  assert diff_sigma_p < 1.e-10, diff_sigma_p
  # Emodel should differ by exactly 1/sqrt(2).
  ratio = flex.abs(r1.e_model) / flex.abs(r2.e_model)
  for v in ratio:
    assert approx_equal(v, 2.0 ** 0.5, eps=1.e-6), v

def exercise_e_model_mean_square_near_one_on_real_fmodel():
  # The defining physical property of a normalised amplitude: mean|E|^2
  # should sit near 1 in every resolution shell, for a real (scaled)
  # fmodel -- not just pass an isolated numerical check on hand-built
  # inputs. Loose tolerance since this is a statistical property on a
  # modest-sized synthetic structure, not an exact identity.
  fmodel = build_fmodel(n_atoms=60, d_min=1.6, seed=6)
  fmnas = llgi_e_bs.f_model_no_aniso_scale(fmodel)
  d_star_sq = fmodel.f_obs().d_star_sq().data()
  epsilons = fmodel.f_obs().epsilons().data().as_double()
  result = llgi_e_bs.build_e_model(fmnas.data(), epsilons, d_star_sq)
  e_model_sq = flex.norm(result.e_model)

  import numpy as np
  order = np.argsort(d_star_sq.as_numpy_array())
  bins = np.array_split(order, 6)
  e_sq_np = e_model_sq.as_numpy_array()
  for b in bins:
    if(len(b) == 0):
      continue
    mean_e_sq = e_sq_np[b].mean()
    assert 0.7 < mean_e_sq < 1.3, mean_e_sq

def exercise_degenerate_resolution_range_does_not_crash():
  # All reflections at (numerically) identical d*^2 -- an edge case
  # _auto_kernel_width / chebyshev fitting must handle without crashing
  # (kernel_normalisation's own auto_kernel heuristic has an explicit
  # fallback loop for this; build_sigma_p relies on that loop, not a
  # separate guard of its own).
  n = 50
  d_star_sq = flex.double([0.1] * n)
  rnd = random.Random(8)
  fmnas = flex.complex_double(
    [complex(rnd.uniform(1.0, 5.0), 0.0) for i in range(n)])
  sigma_p = llgi_e_bs.build_sigma_p(fmnas, d_star_sq, n_nodes=8)
  assert sigma_p.size() == n
  for v in sigma_p:
    assert v > 0

def _synthetic_llgi_inputs(fmodel, seed=10):
  # (dobs, feff, resn) as plain flex.double, on fmodel.f_obs()'s CURRENT
  # index set (i.e. after any outlier removal update_all_scales() did).
  d = synthetic_llgi_data(fmodel, seed=seed)
  return d.dobs.data(), d.feff.data(), d.resn.data()

def exercise_bss_k_sol_b_sol_recovers_known_values():
  # Inject a synthetic k_mask array with a KNOWN (k_sol, b_sol) plus
  # small multiplicative noise into a real fmodel (via
  # fmodel.update(k_mask=...)): bss_k_sol_b_sol should recover both to
  # within a few percent, exercised through the real function against a
  # real fmodel/f_obs, not a hand-rolled duplicate of its fit logic.
  from mmtbx.f_model import ext as f_model_ext
  fmodel = build_fmodel(n_atoms=30, d_min=2.0, seed=14)
  ss = llgi_e_bs.ss_from_f_obs(fmodel.f_obs())
  true_k_sol, true_b_sol = 0.35, 45.0
  k_mask_true = f_model_ext.k_mask(ss, true_k_sol, true_b_sol)
  rnd = random.Random(13)
  n = ss.size()
  noise = flex.double([rnd.gauss(0.0, 0.01) for i in range(n)])
  k_mask_noisy = k_mask_true * (flex.double(n, 1.0) + noise)
  fmodel.update(k_mask=[k_mask_noisy])

  k_sol, b_sol = llgi_e_bs.bss_k_sol_b_sol(fmodel)
  assert approx_equal(k_sol, true_k_sol, eps=0.02), (k_sol, true_k_sol)
  assert approx_equal(b_sol, true_b_sol, eps=1.0), (b_sol, true_b_sol)

def exercise_bss_k_sol_b_sol_falls_back_on_degenerate_input():
  # k_mask ~ 0 everywhere (e.g. a tiny test structure with no real
  # solvent channels) must NOT produce a nonsensical fit (observed
  # directly: an early version gave b_sol ~ -164 on exactly this kind of
  # input) -- bss_k_sol_b_sol should fall back to its defaults.
  fmodel = build_fmodel(n_atoms=15, d_min=2.0, seed=15)
  k_sol, b_sol = llgi_e_bs.bss_k_sol_b_sol(
    fmodel, k_sol_default=0.35, b_sol_default=46.0)
  assert k_sol is not None and b_sol is not None
  assert 0.0 <= k_sol <= 0.6
  assert 0.0 <= b_sol <= 150.0

def exercise_sigmaa_fit_leaves_k_mask_untouched():
  # For both sigmaA forms, estimate_e_sigmaa_for_fmodel must run to
  # completion against a real fmodel, return one sigmaA in (0, 1) per
  # reflection, and leave fmodel's k_mask, exactly as
  # update_all_scales()/bss set it, BIT-IDENTICAL.
  fmodel = build_fmodel(n_atoms=50, d_min=1.75, seed=30)
  dobs, feff, resn = _synthetic_llgi_inputs(fmodel, seed=31)
  k_mask_before = flex.double(fmodel.k_masks()[0])
  for sigmaa_model in ["spline", "d_model"]:
    params = llgi_e_bs.llgi_e_bulk_solvent_params.extract()
    params.sigmaa_model = sigmaa_model
    params.sigmaa_max_iterations = 20
    params.d_model_params.max_iterations = 60
    result = llgi_e_bs.estimate_e_sigmaa_for_fmodel(
      fmodel, dobs=dobs, feff=feff, resn=resn, params=params)
    k_mask_after = flex.double(fmodel.k_masks()[0])
    assert approx_equal(list(k_mask_before), list(k_mask_after), eps=0.0)
    assert result.sigmaa.size() == fmodel.f_obs().size()
    for v in result.sigmaa:
      assert 0.0 < v < 1.0, (sigmaa_model, v)

def exercise_e_sigmaa_curvature_penalty_gradient_finite_difference():
  # End-to-end finite-difference check of e_sigmaa_target_evaluator.
  # compute_functional_and_gradients with the shared spline curvature
  # restraint (mmtbx.refinement.llgi_sigmaa.
  # _spline_curvature_penalty_and_gradient, curvature_weight > 0) mixed
  # into the E-scale LLGI gradient -- mirrors mmtbx.regression.
  # tst_llgi_sigmaa.exercise_curvature_penalty_evaluator_gradient_
  # finite_difference for the F-scale evaluator, since the two
  # evaluators wire the penalty in independently (separate __init__/
  # compute_functional_and_gradients implementations) even though they
  # share the same penalty helper.
  import numpy as np
  fmodel = build_fmodel(n_atoms=40, d_min=1.8, seed=21)
  dobs, feff, resn = _synthetic_llgi_inputs(fmodel, seed=22)
  f_obs = fmodel.f_obs()
  epsilons = f_obs.epsilons().data().as_double()
  d_star_sq = f_obs.d_star_sq().data()
  centric_flags = f_obs.centric_flags().data()

  fmnas = llgi_e_bs.f_model_no_aniso_scale(fmodel).data()
  e_model_cplx = llgi_e_bs.build_e_model(
    fmnas, epsilons, d_star_sq).e_model
  e_model = flex.abs(e_model_cplx)
  e_eff = llgi_e_bs.build_e_eff(feff, resn)

  sigmaa_design = llgi_e_bs._b_spline_design_matrix(
    d_star_sq.as_numpy_array(), 6, 3)
  evaluator = llgi_e_bs.e_sigmaa_target_evaluator(
    e_eff=e_eff, test_selection=fmodel.r_free_flags().data(),
    e_model=e_model, dobs=dobs, centric_flags=centric_flags,
    sigmaa_design=sigmaa_design, n_sigmaa_coeffs=6,
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

def exercise_d_model_sigmaa_is_the_default():
  # sigmaa_model defaults to "d_model" -- pinned here so an accidental
  # flip is caught. The plain extract() (used as the fallback when no
  # params are passed) must carry d_model_params too.
  params = llgi_e_bs.llgi_e_bulk_solvent_params.extract()
  assert params.sigmaa_model == "d_model"
  assert params.d_model_params.include_constant_term is True

def run():
  exercise_f_model_no_aniso_scale_matches_direct_computation()
  exercise_f_model_no_aniso_scale_times_k_aniso_is_f_model()
  exercise_sigma_p_constant_intensity_recovers_constant()
  exercise_sigma_p_is_epsilon_free()
  exercise_e_model_mean_square_near_one_on_real_fmodel()
  exercise_degenerate_resolution_range_does_not_crash()
  exercise_bss_k_sol_b_sol_recovers_known_values()
  exercise_bss_k_sol_b_sol_falls_back_on_degenerate_input()
  exercise_e_sigmaa_curvature_penalty_gradient_finite_difference()
  exercise_sigmaa_fit_leaves_k_mask_untouched()
  exercise_d_model_sigmaa_is_the_default()
  print("OK")

if (__name__ == "__main__"):
  run()
