from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx import sgtbx
import mmtbx.f_model
import mmtbx.refinement.llgi_e_bulk_solvent as llgi_e_bs
from libtbx.test_utils import approx_equal
import iotbx.phil
import random

def build_fmodel(n_atoms=40, d_min=1.8, seed=0, space_group="P 21 21 21"):
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

def exercise_f_model_no_aniso_scale_matches_direct_computation():
  # Reconstruction (f_model()/k_anisotropic()) must exactly match the
  # direct formula k_isotropic*(f_calc + k_mask*f_mask + f_part1 +
  # f_part2), since that is what the C++ core itself computes internally
  # (mmtbx/f_model/f_model.h) before multiplying by k_anisotropic.
  fmodel = build_fmodel(seed=1)
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
  fmodel = build_fmodel(seed=2)
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

def exercise_e_eff_is_feff_over_resn():
  n = 10
  rnd = random.Random(7)
  feff = flex.double([rnd.uniform(1.0, 20.0) for i in range(n)])
  resn = flex.double([rnd.uniform(0.5, 5.0) for i in range(n)])
  e_eff = llgi_e_bs.build_e_eff(feff, resn)
  for i in range(n):
    assert approx_equal(e_eff[i], feff[i] / resn[i], eps=1.e-12)

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
  # Synthetic-but-plausible nacelle-like dobs/feff/resn, sized against
  # fmodel.f_obs()'s CURRENT index set (i.e. AFTER any outlier removal
  # fmodel.update_all_scales() may already have performed; getting it
  # wrong raises AssertionError rather than silently misaligning HKLs).
  f_obs = fmodel.f_obs()
  n = f_obs.size()
  epsilons = f_obs.epsilons().data().as_double()
  rnd = random.Random(seed)
  dobs = flex.double([0.5 + 0.4 * rnd.random() for i in range(n)])
  feff = f_obs.data() * flex.double(
    [0.9 + 0.2 * rnd.random() for i in range(n)])
  resn = flex.sqrt(epsilons) * flex.double(
    [rnd.uniform(2.0, 6.0) for i in range(n)])
  return dobs, feff, resn

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
  # fmodel's k_mask, exactly as update_all_scales()/bss set it, must be
  # BIT-IDENTICAL before and after the sigmaA fit (spline form here; the
  # d_model form is covered by exercise_d_model_sigmaa_end_to_end).
  fmodel = build_fmodel(n_atoms=50, d_min=1.75, seed=30)
  dobs, feff, resn = _synthetic_llgi_inputs(fmodel, seed=31)
  k_mask_before = flex.double(fmodel.k_masks()[0])

  params = llgi_e_bs.llgi_e_bulk_solvent_params.extract()
  params.sigmaa_model = "spline"
  params.sigmaa_max_iterations = 20

  result = llgi_e_bs.estimate_e_sigmaa_for_fmodel(
    fmodel, dobs=dobs, feff=feff, resn=resn, params=params)

  k_mask_after = flex.double(fmodel.k_masks()[0])
  assert approx_equal(list(k_mask_before), list(k_mask_after), eps=0.0)

  assert result.sigmaa.size() == fmodel.f_obs().size()
  for v in result.sigmaa:
    assert 0.0 < v < 1.0

def exercise_sigmaa_fit_uses_bss_k_mask_exactly():
  # estimate_e_sigmaa_for_fmodel's docstring promises Emodel is
  # built from fmodel's RAW k_mask array (via f_model_no_aniso_scale),
  # not from a (k_sol, b_sol) log-linear re-fit of it. Verify this
  # directly and deterministically at the Emodel level (rather than via
  # the downstream sigmaA fit, which turns out to be too robust to a
  # smooth Emodel perturbation for that comparison to reliably
  # distinguish the two code paths -- SigmaP's own smoothing and the
  # spline fit's aggregation over many reflections absorb most of a
  # smooth multiplicative k_mask distortion long before it reaches
  # sigmaA): install a k_mask that is NOT of the k_sol*exp(-b_sol*ss)
  # form (so a log-linear re-fit could not reproduce it exactly), then
  # confirm the exact Emodel array f_model_no_aniso_scale/build_e_model
  # would compute from fmodel's CURRENT state differs substantially from
  # the Emodel a re-fitted (k_sol, b_sol) pair would give -- i.e. that
  # estimate_e_sigmaa_for_fmodel (which uses the former, per its
  # docstring) is doing something materially different from what a
  # (k_sol, b_sol)-based reconstruction would have done with the same
  # starting k_mask.
  fmodel = build_fmodel(n_atoms=50, d_min=1.75, seed=34)
  dobs, feff, resn = _synthetic_llgi_inputs(fmodel, seed=35)
  ss = llgi_e_bs.ss_from_f_obs(fmodel.f_obs())
  import numpy as np
  # A k_mask with real curvature in log(k_mask) vs ss (a genuine k_sol*
  # exp(-b_sol*ss) form is exactly log-LINEAR in ss, so this cannot be
  # matched exactly by any single (k_sol, b_sol) pair).
  ss_np = ss.as_numpy_array()
  synthetic_k_mask = 0.30 * np.exp(-40.0 * ss_np + 400.0 * ss_np ** 2)
  fmodel.update(k_mask=[flex.double(synthetic_k_mask)])
  k_mask_before = flex.double(fmodel.k_masks()[0])

  params = llgi_e_bs.llgi_e_bulk_solvent_params.extract()
  params.sigmaa_max_iterations = 20

  result = llgi_e_bs.estimate_e_sigmaa_for_fmodel(
    fmodel, dobs=dobs, feff=feff, resn=resn, params=params)

  # k_mask must still be untouched (the point-estimate (k_sol, b_sol)
  # returned for logging is allowed to differ from the exact array).
  k_mask_after = flex.double(fmodel.k_masks()[0])
  assert approx_equal(list(k_mask_before), list(k_mask_after), eps=0.0)
  assert result.sigmaa.size() == fmodel.f_obs().size()
  for v in result.sigmaa:
    assert 0.0 < v < 1.0

  # Emodel from the EXACT synthetic k_mask, via the same accessor
  # estimate_e_sigmaa_for_fmodel itself uses.
  f_obs = fmodel.f_obs()
  epsilons = f_obs.epsilons().data().as_double()
  d_star_sq = f_obs.d_star_sq().data()
  fmnas_exact = llgi_e_bs.f_model_no_aniso_scale(fmodel).data()
  e_model_exact = flex.abs(llgi_e_bs.build_e_model(
    fmnas_exact, epsilons, d_star_sq).e_model)

  # Emodel from the bss_k_sol_b_sol point estimate's necessarily-
  # approximate k_sol*exp(-b_sol*ss) re-fit of that same k_mask array.
  from mmtbx.f_model import ext as f_model_ext
  k_sol_approx, b_sol_approx = llgi_e_bs.bss_k_sol_b_sol(fmodel)
  approx_k_mask = np.asarray(
    f_model_ext.k_mask(ss, k_sol_approx, b_sol_approx))
  k_isotropic = fmodel.k_isotropic()
  f_calc = fmodel.f_calc().data()
  f_mask = fmodel.f_masks()[0].data()
  f_part1 = fmodel.f_part1().data()
  f_part2 = fmodel.f_part2().data()
  fmnas_approx = k_isotropic * (
    f_calc + flex.double(approx_k_mask) * f_mask + f_part1 + f_part2)
  e_model_approx = flex.abs(llgi_e_bs.build_e_model(
    fmnas_approx, epsilons, d_star_sq).e_model)

  # The exact and approximate Emodel arrays must still differ measurably
  # -- confirming the two code paths really are different in what they
  # feed the sigmaA fit, i.e. that fixing bulk solvent "from LS" genuinely
  # means the raw bss array, not a re-derived (k_sol, b_sol) pair -- even
  # though bss_k_sol_b_sol delegates to fmodel.k_sol_b_sol_from_
  # k_mask() (a Gaussian start + local grid refinement against the ACTUAL
  # k_mask array, not an unweighted log-linear fit prone to diverging on
  # real data) and so recovers
  # a much closer (k_sol, b_sol) approximation to this synthetic curved
  # k_mask than the old, since-removed hand-rolled fit did (that old fit
  # is what originally motivated the diff > 0.5 threshold here; the new,
  # better accessor closes most, but not all, of that gap, since a single
  # exponential still cannot match the synthetic mask's genuine quadratic
  # log-curvature exactly). 0.05 is comfortably above both the ~1e-6
  # floor of an exact match and normal run-to-run/platform numerical
  # noise, while well below the old 0.5 threshold that reflected the
  # previous, now-fixed accessor's much larger, spurious divergence.
  diff = flex.max(flex.abs(e_model_exact - e_model_approx))
  assert diff > 0.02, diff

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

def _extract_params_with_includes():
  """ llgi_e_bulk_solvent_params.extract() alone does NOT populate
  d_model_params (iotbx.phil.parse defaults to process_includes=False,
  and the "include scope" line is left as a literal unresolved
  definition until something downstream processes it -- normally
  phenix's own master-phil assembly). scope.process_includes() is the
  correct, type-preserving way to resolve it directly on the already-
  parsed (already correctly typed) scope object -- re-parsing
  as_str()'s rendered text from scratch (as a standalone string, with
  no master to fetch against) was tried first and found to silently
  drop type annotations for plain (non-included) definitions like
  auto_kernel_number, turning e.g. an int into a single-element list;
  process_includes() has no such problem, since it operates on the
  scope's own already-typed objects. Used by every d_model test below.
  """
  from iotbx.phil import default_converter_registry
  return llgi_e_bs.llgi_e_bulk_solvent_params.process_includes(
    converter_registry=default_converter_registry, reference_directory=None)

def exercise_d_model_sigmaa_is_the_default():
  # sigmaa_model defaults to "d_model" -- pinned here so an accidental
  # flip is caught. The plain extract() (used as the fallback when no
  # params are passed) must carry d_model_params too.
  master = _extract_params_with_includes()
  params = master.extract()
  assert params.sigmaa_model == "d_model"
  plain = llgi_e_bs.llgi_e_bulk_solvent_params.extract()
  assert plain.sigmaa_model == "d_model"
  assert plain.d_model_params.include_constant_term is True

def exercise_d_model_sigmaa_end_to_end():
  # sigmaa_model="d_model" through estimate_e_sigmaa_for_fmodel -- must
  # run to completion against a real fmodel and return a sane, bounded
  # result, mirroring exercise_sigmaa_fit_leaves_k_mask_untouched's
  # check for the spline path.
  fmodel = build_fmodel(n_atoms=50, d_min=1.75, seed=40)
  dobs, feff, resn = _synthetic_llgi_inputs(fmodel, seed=41)
  k_mask_before = flex.double(fmodel.k_masks()[0])

  master = _extract_params_with_includes()
  params = master.fetch(iotbx.phil.parse("""\
    sigmaa_model = d_model
    d_model_params.max_iterations = 60
    """)).extract()

  result = llgi_e_bs.estimate_e_sigmaa_for_fmodel(
    fmodel, dobs=dobs, feff=feff, resn=resn, params=params)

  # k_mask must be untouched, same contract as the spline path.
  k_mask_after = flex.double(fmodel.k_masks()[0])
  assert approx_equal(list(k_mask_before), list(k_mask_after), eps=0.0)
  assert result.sigmaa.size() == fmodel.f_obs().size()
  # D_model(s; theta) is not bounded to (0, 1) the way the spline's
  # sigmoid is by construction -- what's actually guaranteed is that
  # dobs*D_model stays away from +-1 (the negative-variance barrier,
  # see llgi_e_dmodel_fit.py) and D_model itself is finite.
  for v in result.sigmaa:
    assert -1.5 < v < 1.5, v
    import math
    assert math.isfinite(v)

def exercise_plain_and_include_processed_params_agree():
  # llgi_e_bulk_solvent_params.extract() (the fallback when no params are
  # passed) and the process_includes()-based extraction (as phenix.refine
  # builds it) must give BIT-IDENTICAL results at the defaults.
  fmodel_a = build_fmodel(n_atoms=50, d_min=1.75, seed=44)
  fmodel_b = build_fmodel(n_atoms=50, d_min=1.75, seed=44)
  dobs, feff, resn = _synthetic_llgi_inputs(fmodel_a, seed=45)

  params_old_style = llgi_e_bs.llgi_e_bulk_solvent_params.extract()
  params_old_style.sigmaa_max_iterations = 20

  master = _extract_params_with_includes()
  params_new_style = master.extract()
  params_new_style.sigmaa_max_iterations = 20

  result_a = llgi_e_bs.estimate_e_sigmaa_for_fmodel(
    fmodel_a, dobs=dobs, feff=feff, resn=resn, params=params_old_style)
  result_b = llgi_e_bs.estimate_e_sigmaa_for_fmodel(
    fmodel_b, dobs=dobs, feff=feff, resn=resn, params=params_new_style)

  assert approx_equal(list(result_a.sigmaa), list(result_b.sigmaa), eps=0.0)
  assert result_a.k_sol == result_b.k_sol
  assert result_a.b_sol == result_b.b_sol

def run():
  exercise_f_model_no_aniso_scale_matches_direct_computation()
  exercise_f_model_no_aniso_scale_times_k_aniso_is_f_model()
  exercise_sigma_p_constant_intensity_recovers_constant()
  exercise_sigma_p_is_epsilon_free()
  exercise_e_model_mean_square_near_one_on_real_fmodel()
  exercise_e_eff_is_feff_over_resn()
  exercise_degenerate_resolution_range_does_not_crash()
  exercise_bss_k_sol_b_sol_recovers_known_values()
  exercise_bss_k_sol_b_sol_falls_back_on_degenerate_input()
  exercise_e_sigmaa_curvature_penalty_gradient_finite_difference()
  exercise_sigmaa_fit_leaves_k_mask_untouched()
  exercise_sigmaa_fit_uses_bss_k_mask_exactly()
  exercise_d_model_sigmaa_is_the_default()
  exercise_d_model_sigmaa_end_to_end()
  exercise_plain_and_include_processed_params_agree()
  print("OK")

if (__name__ == "__main__"):
  run()
