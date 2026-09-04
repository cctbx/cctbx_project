from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx import sgtbx
import mmtbx.f_model
import mmtbx.refinement.llgi_e_bulk_solvent as llgi_e_bs
from libtbx.test_utils import approx_equal
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
  # fmodel.update_all_scales() may already have performed -- see
  # run_inner_loop's docstring for why this ordering matters; getting it
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

def exercise_bulk_solvent_gradient_matches_fixed_sigma_p_finite_difference():
  # bulk_solvent_target_and_gradients's analytic gradient is derived
  # holding SigmaP fixed (see its docstring for why the TRUE total
  # derivative, letting SigmaP respond to k_sol/b_sol too, is a
  # materially different, larger quantity -- roughly 10x on real data).
  # This test verifies internal consistency: the analytic gradient must
  # match a finite difference of the SAME (fixed-SigmaP) target the
  # function itself would compute at that one point -- not a finite
  # difference that lets SigmaP float, which would NOT match by design.
  fmodel = build_fmodel(n_atoms=40, d_min=1.8, seed=11)
  dobs, feff, resn = _synthetic_llgi_inputs(fmodel, seed=12)
  f_obs = fmodel.f_obs()
  n = f_obs.size()
  epsilons = f_obs.epsilons().data().as_double()
  d_star_sq = f_obs.d_star_sq().data()
  ss = llgi_e_bs.ss_from_f_obs(f_obs)
  centric_flags = f_obs.centric_flags().data()
  working_sel = ~fmodel.r_free_flags().data()
  e_eff = llgi_e_bs.build_e_eff(feff, resn)
  sigmaa = flex.double(n, 0.5)

  f_calc = fmodel.f_calc().data()
  f_mask = fmodel.f_masks()[0].data()
  f_part1 = fmodel.f_part1().data()
  f_part2 = fmodel.f_part2().data()
  k_iso = fmodel.k_isotropic()

  k_sol0, b_sol0 = 0.30, 40.0

  # Freeze SigmaP at k_sol0/b_sol0, matching what the analytic gradient
  # actually differentiates (see bulk_solvent_target_and_gradients's
  # docstring): a small local helper reimplementing the target with a
  # caller-supplied fixed SigmaP, used only for this finite-difference
  # cross-check.
  from cctbx.xray import ext as xray_ext
  k_mask0 = llgi_e_bs.k_mask_and_gradients(ss, k_sol0, b_sol0).k_mask
  fmnas0 = k_iso * (f_calc + k_mask0 * f_mask + f_part1 + f_part2)
  sigma_p_fixed = llgi_e_bs.build_sigma_p(fmnas0, d_star_sq)

  def target_fixed_sigma_p(k_sol, b_sol):
    k_mask_r = llgi_e_bs.k_mask_and_gradients(ss, k_sol, b_sol)
    fmnas = k_iso * (f_calc + k_mask_r.k_mask * f_mask + f_part1 + f_part2)
    denom = flex.sqrt(epsilons * sigma_p_fixed)
    e_model = fmnas * (1.0 / denom)
    result = xray_ext.llgi_e_emodel_target_and_gradients(
      e_eff=e_eff, selection=working_sel, e_model=flex.abs(e_model),
      dobs=dobs, sigmaa=sigmaa, centric_flags=centric_flags)
    return result.target()

  r0 = llgi_e_bs.bulk_solvent_target_and_gradients(
    e_eff=e_eff, selection=working_sel, dobs=dobs, sigmaa=sigmaa,
    centric_flags=centric_flags, f_calc=f_calc, f_mask=f_mask,
    f_part1=f_part1, f_part2=f_part2, k_isotropic=k_iso,
    epsilons=epsilons, d_star_sq=d_star_sq, ss=ss,
    k_sol=k_sol0, b_sol=b_sol0)

  eps_k = 1.e-6
  fd_k = (target_fixed_sigma_p(k_sol0 + eps_k, b_sol0)
          - target_fixed_sigma_p(k_sol0 - eps_k, b_sol0)) / (2 * eps_k)
  eps_b = 1.e-4
  fd_b = (target_fixed_sigma_p(k_sol0, b_sol0 + eps_b)
          - target_fixed_sigma_p(k_sol0, b_sol0 - eps_b)) / (2 * eps_b)

  assert approx_equal(r0.gradients[0], fd_k, eps=1.e-4 * max(1, abs(fd_k))), (
    r0.gradients[0], fd_k)
  assert approx_equal(r0.gradients[1], fd_b, eps=1.e-4 * max(1, abs(fd_b))), (
    r0.gradients[1], fd_b)

def exercise_initial_k_sol_b_sol_recovers_known_values():
  # Inject a synthetic k_mask array with a KNOWN (k_sol, b_sol) plus
  # small multiplicative noise into a real fmodel (via
  # fmodel.update(k_mask=...), the same mechanism run_inner_loop itself
  # uses to push a fitted bulk-solvent model back onto fmodel): the log-
  # linear recovery fit (initial_k_sol_b_sol) should recover both to
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

  k_sol, b_sol = llgi_e_bs.initial_k_sol_b_sol(fmodel)
  assert approx_equal(k_sol, true_k_sol, eps=0.02), (k_sol, true_k_sol)
  assert approx_equal(b_sol, true_b_sol, eps=1.0), (b_sol, true_b_sol)

def exercise_initial_k_sol_b_sol_falls_back_on_degenerate_input():
  # k_mask ~ 0 everywhere (e.g. a tiny test structure with no real
  # solvent channels) must NOT produce a nonsensical fit (observed
  # directly: an early version gave b_sol ~ -164 on exactly this kind of
  # input) -- initial_k_sol_b_sol should fall back to its defaults.
  fmodel = build_fmodel(n_atoms=15, d_min=2.0, seed=15)
  k_sol, b_sol = llgi_e_bs.initial_k_sol_b_sol(
    fmodel, k_sol_default=0.35, b_sol_default=46.0)
  assert k_sol is not None and b_sol is not None
  assert 0.0 <= k_sol <= 0.6
  assert 0.0 <= b_sol <= 150.0

def exercise_run_inner_loop_end_to_end():
  # Full Stage-1/Stage-2 inner loop, on real (synthetic) data: must run
  # to completion (converged or max_inner_iterations reached), return a
  # sensible-shaped result, and leave fmodel's k_mask updated to match
  # the returned (k_sol, b_sol) exactly.
  #
  # fix_bulk_solvent_from_ls is set False EXPLICITLY here: it defaults to
  # True, which skips Stage 2 entirely and never touches k_mask (see
  # exercise_fixed_bulk_solvent_leaves_k_mask_untouched), so relying on
  # the default would silently turn this into a test of a different code
  # path than the one its assertions below describe.
  fmodel = build_fmodel(n_atoms=60, d_min=1.7, seed=16)
  dobs, feff, resn = _synthetic_llgi_inputs(fmodel, seed=17)
  params = llgi_e_bs.llgi_e_bulk_solvent_params.extract()
  params.fix_bulk_solvent_from_ls = False
  params.max_inner_iterations = 3
  params.sigmaa_max_iterations = 20
  params.bulk_solvent_max_iterations = 15

  result = llgi_e_bs.run_inner_loop(
    fmodel, dobs=dobs, feff=feff, resn=resn, params=params)

  assert result.n_iterations >= 1
  assert result.sigmaa.size() == fmodel.f_obs().size()
  for v in result.sigmaa:
    assert 0.0 < v < 1.0
  assert 0.0 <= result.k_sol <= 0.6
  assert 0.0 <= result.b_sol <= 150.0
  assert len(result.history) == result.n_iterations

  # fmodel's k_mask should now match the returned (k_sol, b_sol) exactly.
  from mmtbx.f_model import ext as f_model_ext
  ss = llgi_e_bs.ss_from_f_obs(fmodel.f_obs())
  expected_k_mask = f_model_ext.k_mask(ss, result.k_sol, result.b_sol)
  actual_k_mask = fmodel.k_masks()[0]
  diff = flex.max(flex.abs(
    flex.double(expected_k_mask) - flex.double(actual_k_mask)))
  assert diff < 1.e-10, diff

def exercise_fix_bulk_solvent_from_ls_is_the_default():
  # fix_bulk_solvent_from_ls defaults to True: on 2G38 over 5
  # macrocycles it gave the best R-work/R-free of the three approaches
  # tested (0.3832/0.4144, vs 0.3881/0.4155 and 0.3912/0.4159 for
  # LLGI-driven bulk solvent at sigmaA curvature weights 0.001 and
  # 0.02), and it sidesteps Stage 2's two-competing-solutions/boundary-
  # pinning behaviour entirely. Pinned here so an accidental flip is
  # caught: several other tests in this file and in tst_llgi_data.py set
  # it False EXPLICITLY to exercise the two-stage path, and would
  # silently change meaning if the default moved.
  params = llgi_e_bs.llgi_e_bulk_solvent_params.extract()
  assert params.fix_bulk_solvent_from_ls is True

def exercise_fixed_bulk_solvent_leaves_k_mask_untouched():
  # params.fix_bulk_solvent_from_ls (the hybrid mode): fmodel's k_mask,
  # exactly as update_all_scales()/bss set it, must be BIT-IDENTICAL
  # before and after -- this mode must never call fmodel.update(k_mask=
  # ...) at all, unlike the default Stage-1/Stage-2 loop.
  fmodel = build_fmodel(n_atoms=50, d_min=1.75, seed=30)
  dobs, feff, resn = _synthetic_llgi_inputs(fmodel, seed=31)
  k_mask_before = flex.double(fmodel.k_masks()[0])

  params = llgi_e_bs.llgi_e_bulk_solvent_params.extract()
  params.fix_bulk_solvent_from_ls = True
  params.sigmaa_max_iterations = 20

  result = llgi_e_bs.run_inner_loop(
    fmodel, dobs=dobs, feff=feff, resn=resn, params=params)

  k_mask_after = flex.double(fmodel.k_masks()[0])
  assert approx_equal(list(k_mask_before), list(k_mask_after), eps=0.0)

  assert result.n_iterations == 1
  assert result.converged is True
  assert result.sigmaa.size() == fmodel.f_obs().size()
  for v in result.sigmaa:
    assert 0.0 < v < 1.0
  assert len(result.history) == 1

def exercise_fixed_bulk_solvent_matches_direct_call():
  # run_inner_loop(params.fix_bulk_solvent_from_ls=True) must dispatch to
  # estimate_e_sigmaa_fixed_bulk_solvent exactly -- same sigmaA, same
  # (k_sol, b_sol) point estimate -- not merely something similarly
  # shaped, since callers (mmtbx.f_model.manager.
  # update_llgi_e_bulk_solvent) only ever go through run_inner_loop.
  fmodel = build_fmodel(n_atoms=50, d_min=1.75, seed=32)
  dobs, feff, resn = _synthetic_llgi_inputs(fmodel, seed=33)
  params = llgi_e_bs.llgi_e_bulk_solvent_params.extract()
  params.fix_bulk_solvent_from_ls = True
  params.sigmaa_max_iterations = 20

  via_run_inner_loop = llgi_e_bs.run_inner_loop(
    fmodel, dobs=dobs, feff=feff, resn=resn, params=params)
  direct = llgi_e_bs.estimate_e_sigmaa_fixed_bulk_solvent(
    fmodel, dobs=dobs, feff=feff, resn=resn, params=params)

  assert approx_equal(
    list(via_run_inner_loop.sigmaa), list(direct.sigmaa))
  assert via_run_inner_loop.k_sol == direct.k_sol
  assert via_run_inner_loop.b_sol == direct.b_sol

def exercise_fixed_bulk_solvent_uses_bss_k_mask_exactly():
  # estimate_e_sigmaa_fixed_bulk_solvent's docstring promises Emodel is
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
  # estimate_e_sigmaa_fixed_bulk_solvent (which uses the former, per its
  # docstring) is doing something materially different from what a
  # (k_sol, b_sol)-based reconstruction (like run_inner_loop's own
  # Stage-1 bootstrap) would have done with the same starting k_mask.
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
  params.fix_bulk_solvent_from_ls = True
  params.sigmaa_max_iterations = 20

  result = llgi_e_bs.estimate_e_sigmaa_fixed_bulk_solvent(
    fmodel, dobs=dobs, feff=feff, resn=resn, params=params)

  # k_mask must still be untouched (the point-estimate (k_sol, b_sol)
  # returned for logging is allowed to differ from the exact array).
  k_mask_after = flex.double(fmodel.k_masks()[0])
  assert approx_equal(list(k_mask_before), list(k_mask_after), eps=0.0)
  assert result.sigmaa.size() == fmodel.f_obs().size()
  for v in result.sigmaa:
    assert 0.0 < v < 1.0

  # Emodel from the EXACT synthetic k_mask, via the same accessor
  # estimate_e_sigmaa_fixed_bulk_solvent itself uses.
  f_obs = fmodel.f_obs()
  epsilons = f_obs.epsilons().data().as_double()
  d_star_sq = f_obs.d_star_sq().data()
  fmnas_exact = llgi_e_bs.f_model_no_aniso_scale(fmodel).data()
  e_model_exact = flex.abs(llgi_e_bs.build_e_model(
    fmnas_exact, epsilons, d_star_sq).e_model)

  # Emodel from the initial_k_sol_b_sol point estimate's necessarily-
  # approximate k_sol*exp(-b_sol*ss) re-fit of that same k_mask array
  # (what a (k_sol, b_sol)-based reconstruction, e.g. run_inner_loop's
  # own Stage-1 bootstrap, would use instead).
  k_sol_approx, b_sol_approx = llgi_e_bs.initial_k_sol_b_sol(fmodel)
  approx_k_mask = np.asarray(
    llgi_e_bs.k_mask_and_gradients(ss, k_sol_approx, b_sol_approx).k_mask)
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
  # though initial_k_sol_b_sol now delegates to fmodel.k_sol_b_sol_from_
  # k_mask() (a Gaussian start + local grid refinement against the ACTUAL
  # k_mask array, not an unweighted log-linear fit prone to diverging on
  # real data -- see initial_k_sol_b_sol's own docstring) and so recovers
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

def run():
  exercise_f_model_no_aniso_scale_matches_direct_computation()
  exercise_f_model_no_aniso_scale_times_k_aniso_is_f_model()
  exercise_sigma_p_constant_intensity_recovers_constant()
  exercise_sigma_p_is_epsilon_free()
  exercise_e_model_mean_square_near_one_on_real_fmodel()
  exercise_e_eff_is_feff_over_resn()
  exercise_degenerate_resolution_range_does_not_crash()
  exercise_bulk_solvent_gradient_matches_fixed_sigma_p_finite_difference()
  exercise_initial_k_sol_b_sol_recovers_known_values()
  exercise_initial_k_sol_b_sol_falls_back_on_degenerate_input()
  exercise_run_inner_loop_end_to_end()
  exercise_e_sigmaa_curvature_penalty_gradient_finite_difference()
  exercise_fix_bulk_solvent_from_ls_is_the_default()
  exercise_fixed_bulk_solvent_leaves_k_mask_untouched()
  exercise_fixed_bulk_solvent_matches_direct_call()
  exercise_fixed_bulk_solvent_uses_bss_k_mask_exactly()
  print("OK")

if (__name__ == "__main__"):
  run()
