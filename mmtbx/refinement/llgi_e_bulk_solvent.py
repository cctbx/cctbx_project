from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.xray import ext as xray_ext
from mmtbx import scaling
from mmtbx.f_model import ext as f_model_ext
from mmtbx.refinement.llgi_sigmaa import _b_spline_design_matrix
from mmtbx.refinement.llgi_sigmaa import _spline_curvature_penalty_and_gradient
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit
import scitbx.lbfgs
from libtbx import group_args
import iotbx.phil

llgi_e_bulk_solvent_params = iotbx.phil.parse("""\
  n_sigmap_nodes = 15
    .type = int
    .expert_level = 3
    .help = "Number of Chebyshev nodes (and fit terms) used for the " \
            "SigmaP(d*^2) smoothed intensity trend -- see " \
            "mmtbx.refinement.llgi_e_bulk_solvent.build_sigma_p. Matches " \
            "mmtbx.scaling.absolute_scaling.kernel_normalisation's " \
            "n_bins/n_term defaults (23/13); kept smaller and equal here " \
            "since this is recomputed every inner-loop iteration (see " \
            "doc/llgi_target_design.md's E-scale design note sec. 4/7), " \
            "not once per dataset."
  auto_kernel_number = 50
    .type = int
    .expert_level = 3
    .help = "Number of sorted reflections spanned by the auto-tuned " \
            "kernel width, matching kernel_normalisation's own " \
            "number_of_sorted_reflections_for_auto_kernel default."
  n_sigmaa_coeffs = 8
    .type = int
    .short_caption = Number of sigmaA spline coefficients (E-scale)
    .help = "Number of B-spline coefficients for the E-scale " \
            "sigmaA(resolution) curve fit against R-free data each " \
            "inner-loop iteration (design note sec. 5). Same knot-" \
            "placement convention as the F-scale llgi_sigmaa module."
  spline_degree = 3
    .type = int
    .expert_level = 3
    .help = "Degree of the clamped B-spline basis used for sigmaA(d*^2)."
  sigmaa_max_iterations = 100
    .type = int
    .expert_level = 3
    .help = "Max LBFGS iterations for each inner-loop sigmaA refit."
  sigmaa_curvature_weight = 0.02
    .type = float
    .short_caption = SigmaA spline curvature restraint weight (E-scale)
    .help = "Weight of a light restraint on the second difference of the "\
            "E-scale sigmaA B-spline's raw (pre-sigmoid) coefficients, "\
            "penalising curvature of the fitted curve. Same restraint, "\
            "same rationale, as mmtbx.refinement.llgi_sigmaa."\
            "llgi_sigmaa_scatfrac_params.sigmaa_curvature_weight -- see "\
            "mmtbx.refinement.llgi_sigmaa."\
            "_spline_curvature_penalty_and_gradient. 0 disables it."
  k_sol_min = 0.0
    .type = float
    .expert_level = 3
  k_sol_max = 0.6
    .type = float
    .expert_level = 3
  b_sol_min = 0.0
    .type = float
    .expert_level = 3
  b_sol_max = 150.0
    .type = float
    .expert_level = 3
  bulk_solvent_max_iterations = 50
    .type = int
    .expert_level = 3
    .help = "Max LBFGS iterations for each inner-loop bulk-solvent refit."
  max_inner_iterations = 10
    .type = int
    .help = "Max number of Stage-1/Stage-2 (sigmaA / bulk-solvent) " \
            "alternations within one macrocycle (design note sec. 7). " \
            "The loop may stop earlier on convergence (see " \
            "convergence_tolerance)."
  convergence_tolerance = 1.e-4
    .type = float
    .expert_level = 3
    .help = "Inner loop stops early once the working-set LLG (E-scale, " \
            "Stage 2's own objective) changes by less than this between " \
            "consecutive iterations."
  fix_bulk_solvent_from_ls = True
    .type = bool
    .short_caption = Fix bulk solvent at the LS (bss) value, fit only sigmaA
    .help = "DEFAULT (True): keep bss's own least-squares bulk-solvent " \
            "result exactly as it is and only run the LLGI " \
            "sigmaA(resolution) fit (Stage 1) against it, instead of " \
            "letting the E-scale LLGI target drive (k_sol, b_sol) at all " \
            "(Stage 2). run_inner_loop then calls " \
            "estimate_e_sigmaa_fixed_bulk_solvent instead of running its " \
            "own Stage-1/Stage-2 loop; k_sol_min/max, b_sol_min/max, " \
            "bulk_solvent_max_iterations, max_inner_iterations and " \
            "convergence_tolerance are all unused (Stage 2 never runs), " \
            "and fmodel's k_mask is left exactly as bss set it, never " \
            "overwritten. Set False to restore the original two-stage " \
            "design (LLGI drives bulk solvent too). True is the default " \
            "because Stage 2 was repeatedly observed driving (k_sol, " \
            "b_sol) to one of two competing solutions with one or both " \
            "pinned at a phil boundary, from the very first macrocycle; " \
            "a sigmaA-spline curvature restraint (see " \
            "sigmaa_curvature_weight) helps at macrocycle 1 but stops " \
            "mattering by macrocycle 2, and an attempt to derive its " \
            "weight from the LLGI target's own local curvature failed by " \
            "2-3 orders of magnitude. On 2G38 over 5 macrocycles this " \
            "mode gave the best result of the three tested (R-work/" \
            "R-free 0.3832/0.4144, vs 0.3881/0.4155 and 0.3912/0.4159 " \
            "for LLGI-driven bulk solvent at curvature weights 0.001 and " \
            "0.02) -- see doc/llgi_target_design.md sec. 6's LS-fixed " \
            "bulk solvent, LLGI-only sigmaA note. Note this gives up (for " \
            "now) the original idea that LLGI might improve on the LS " \
            "bulk-solvent fit; that was not tested on a case where LS " \
            "itself struggles, so Stage 2 is kept, not deleted."
""")

def _auto_kernel_width(d_star_sq, number=50):
  """ Reproduce mmtbx.scaling.absolute_scaling.kernel_normalisation's
  auto_kernel=True width heuristic exactly (same algorithm, including its
  fallback loop for the widely-spread-but-locally-degenerate case), so
  SigmaP uses the identical auto-tuning the design note's decision (reuse
  kernel_normalisation's auto-tuning, no separate bandwidth parameter)
  calls for -- without going through the kernel_normalisation class
  itself, which divides out epsilon internally (see build_sigma_p).

  Unlike the original (called once per dataset against real experimental
  data, where every reflection sharing exactly the same d*^2 is
  vanishingly unlikely), this is called every inner-loop iteration (see
  the design note sec. 7), so it adds one guard the original does not
  have: if d_star_sq is uniformly (or near-uniformly, within floating-
  point noise) constant across ALL reflections, the original's fallback
  loop can never find a positive width and would raise AssertionError
  (reproduced and caught by tst_llgi_e_bulk_solvent.py's
  exercise_degenerate_resolution_range_does_not_crash). In that case
  there is no meaningful "resolution trend" to smooth over anyway, so
  fall back to a nominal small positive width instead of crashing.
  """
  import numpy as np
  d_star_sq_np = np.asarray(d_star_sq, dtype=float)
  assert d_star_sq_np.size > 1
  d_star_sq_low = float(d_star_sq_np.min())
  d_star_sq_high = float(d_star_sq_np.max())
  if(d_star_sq_high - d_star_sq_low < 1.e-12):
    # All reflections at (numerically) the same resolution: no trend to
    # smooth over. Return a nominal width so the caller's Chebyshev fit
    # still runs (every node ends up seeing the same local average).
    return 1.0
  sort_permut = np.argsort(d_star_sq_np)
  number = min(number, d_star_sq_np.size - 1)
  kernel_width = d_star_sq_np[sort_permut[number]] - d_star_sq_low
  if(kernel_width == 0):
    original_number = number
    while number < d_star_sq_np.size / 2:
      number += original_number
      number = min(number, d_star_sq_np.size - 1)
      kernel_width = d_star_sq_np[sort_permut[number]] - d_star_sq_low
      if(kernel_width > 0):
        break
  if(kernel_width <= 0):
    # Fallback loop exhausted without finding a positive width (e.g. a
    # large majority of reflections share the same d*^2, with only a
    # small tail spread out) -- use the full data range as a last
    # resort rather than crash.
    kernel_width = d_star_sq_high - d_star_sq_low
  assert kernel_width > 0
  return kernel_width

def build_sigma_p(f_model_no_aniso_scale, d_star_sq, n_nodes=15,
      auto_kernel_number=50, d_star_sq_eval=None):
  """ SigmaP(d*^2): a smoothed, EPSILON-FREE resolution trend of
  |f_model_no_aniso_scale|^2, per doc/llgi_target_design.md's E-scale
  design note sec. 4. Deliberately calls the low-level mmtbx.scaling.
  kernel_normalisation extension function directly, with epsilon forced
  to all-ones, rather than going through the kernel_normalisation Python
  class (mmtbx/scaling/absolute_scaling.py): that class divides by the
  per-reflection epsilon internally, before smoothing (absolute_scaling.h:
  norma_I_array[ii] += I_hkl[jj]*result/epsilon_hkl[jj]) -- exactly the
  step this function must NOT take, so the caller (build_e_model) can
  reapply EPS explicitly, per reflection, afterward. An earlier draft of
  this design fed f_model_no_aniso_scale through the Python class
  directly, silently losing that separation; corrected once traced back
  to absolute_scaling.h's source.

  Otherwise follows kernel_normalisation's own recipe exactly: Gaussian-
  kernel-weighted local average of |f_model_no_aniso_scale|^2 (epsilon
  suppressed) at Chebyshev nodes in d*^2, log-space Chebyshev polynomial
  fit through those node values, evaluated back at every reflection's own
  d*^2 -- genuinely smooth (a fitted polynomial), not a step function.

  Kernel width is auto-tuned via _auto_kernel_width, reproducing
  kernel_normalisation's auto_kernel=True heuristic exactly (see the
  design note sec. 9: reuse the existing auto-tuning as-is, no separate
  bandwidth parameter).

  Returns a flex.double, one SigmaP value per input reflection (evaluated
  at that reflection's own d_star_sq, or at d_star_sq_eval instead if
  given -- see below). Recomputed once per inner-loop iteration (sec. 7
  of the design note), since f_model_no_aniso_scale changes with the
  bulk-solvent model.

  d_star_sq_eval: if given, evaluate the FITTED curve (fit against
  f_model_no_aniso_scale/d_star_sq as usual) at this DIFFERENT set of
  d_star_sq values instead of the ones the curve was fit against -- for
  a genuinely disjoint reflection set with no f_model_no_aniso_scale of
  its own (e.g. mmtbx.map_tools.model_missing_reflections_llgi's missing
  -reflection map-coefficient fill, which has no observed data at all to
  build a fresh SigmaP fit from). Clamped to [d_star_sq.min(),
  d_star_sq.max()] before evaluating -- unlike the B-spline sigmaA fit's
  own scipy extrapolate=False (silently zero outside the fit range),
  chebyshev_polynome IS a genuine polynomial and will extrapolate
  without complaint if asked to, which for a resolution trend fit this
  way could blow up arbitrarily outside the range it was actually
  constrained by data; clamping to the nearest endpoint's fitted value
  avoids that failure mode, matching e_sigmaa_target_evaluator.
  evaluate_at's own clamping convention.
  """
  d_star_sq_np = d_star_sq.as_numpy_array()
  d_star_sq_low = float(d_star_sq_np.min())
  d_star_sq_high = float(d_star_sq_np.max())
  intensity = flex.norm(f_model_no_aniso_scale)  # |fmnas|^2
  epsilon_ones = flex.double(intensity.size(), 1.0)
  kernel_width = _auto_kernel_width(
    d_star_sq_np, number=auto_kernel_number)
  nodes = chebyshev_lsq_fit.chebyshev_nodes(
    n=n_nodes, low=d_star_sq_low, high=d_star_sq_high, include_limits=True)
  mean_intensity_at_nodes = scaling.kernel_normalisation(
    d_star_sq_hkl=d_star_sq,
    I_hkl=intensity,
    epsilon=epsilon_ones,
    d_star_sq_array=nodes,
    kernel_width=kernel_width)
  # Guard against non-positive smoothed values (e.g. a node outside the
  # data's support, or numerical noise at the tails) before taking log,
  # matching kernel_normalisation's own eps=1e-16-style floor.
  floor = 1.e-16
  mean_intensity_at_nodes = flex.double([
    max(v, floor) for v in mean_intensity_at_nodes])
  log_mean_at_nodes = flex.log(mean_intensity_at_nodes)
  fit = chebyshev_lsq_fit.chebyshev_lsq_fit(
    n_nodes, nodes, log_mean_at_nodes)
  poly = chebyshev_polynome(
    n_nodes, d_star_sq_low, d_star_sq_high, fit.coefs)
  if(d_star_sq_eval is None):
    eval_at = d_star_sq
  else:
    import numpy as np
    eval_np = np.clip(
      d_star_sq_eval.as_numpy_array(), d_star_sq_low, d_star_sq_high)
    eval_at = flex.double(eval_np.tolist())
  fitted_log = poly.f(eval_at)
  return flex.exp(fitted_log)

def build_e_model(f_model_no_aniso_scale, epsilons, d_star_sq,
      n_sigmap_nodes=15, auto_kernel_number=50):
  """ Emodel = f_model_no_aniso_scale / sqrt(EPS * SigmaP), per the design
  note sec. 4 (as corrected: EPS re-introduced explicitly, per reflection,
  after SigmaP's smoothing step, which is deliberately epsilon-free -- see
  build_sigma_p).

  epsilons: per-reflection space-group multiplicity/epsilon factor (e.g.
  f_obs.epsilons().data().as_double()).

  Returns a group_args with .e_model (flex.complex_double, same order as
  the input) and .sigma_p (flex.double, for diagnostics/logging -- e.g.
  comparing against the F-scale ScatFrac(resolution) curve).
  """
  n = f_model_no_aniso_scale.size()
  assert epsilons.size() == n
  assert d_star_sq.size() == n
  sigma_p = build_sigma_p(
    f_model_no_aniso_scale, d_star_sq,
    n_nodes=n_sigmap_nodes, auto_kernel_number=auto_kernel_number)
  denom = flex.sqrt(epsilons * sigma_p)
  inv_denom = 1.0 / denom
  e_model = f_model_no_aniso_scale * inv_denom
  return group_args(e_model=e_model, sigma_p=sigma_p)

def build_e_eff(feff, resn):
  """ Eeff = Feff / RESN, per the design note sec. 4: RESN is nacelle's
  own "Root-EpsilonSigmaN" normaliser, already epsilon- and Wilson-trend-
  corrected via nacelle's own Bayesian modelling of the expected
  intensity over reciprocal space (including whatever anisotropy/tNCS
  treatment nacelle applies on the experimental side). No further
  smoothing/normalisation step is applied here -- re-deriving it via a
  second, independent kernel-smoothing pass would be redundant at best
  and inconsistent with nacelle's own normalisation at worst (an earlier
  draft of this design proposed exactly that second pass; corrected once
  RESN's actual meaning -- "Root-EpsilonSigmaN" -- was clarified).

  Computed once per macrocycle (not once per inner-loop iteration -- the
  experimental data does not change within the inner loop; see the
  design note sec. 4/7), directly from the nacelle FEFF/RESN columns
  already loaded via phenix.refinement.llgi_data.get_llgi_data.

  Returns a flex.double, one Eeff value per input reflection.
  """
  assert resn.size() == feff.size()
  return feff / resn

def _sigmoid(z, lower=0.01, upper=0.99):
  """ Bounded sigmoid reparameterisation, identical convention to
  mmtbx.refinement.llgi_sigmaa._sigmoid (same lower/upper bounds), so the
  E-scale sigmaA fit below matches the F-scale one's behaviour exactly.
  Returns (value, d(value)/d(z)), both numpy arrays.
  """
  import numpy as np
  z = np.asarray(z, dtype=float)
  exp_neg_z = np.exp(-z)
  value = lower + (upper - lower) / (1.0 + exp_neg_z)
  dvalue_dz = (upper - lower) * exp_neg_z / (1.0 + exp_neg_z) ** 2
  return value, dvalue_dz

def _bounded(z, lower, upper):
  """ Bounded sigmoid reparameterisation with caller-supplied bounds, for
  k_sol/b_sol (whose natural ranges are very different from sigmaA's
  (0,1) and from each other -- see llgi_e_bulk_solvent_params.k_sol_min/
  max, b_sol_min/max). Same functional form as _sigmoid, generalised.
  Returns (value, d(value)/d(z)).
  """
  import numpy as np
  z = np.asarray(z, dtype=float)
  exp_neg_z = np.exp(-z)
  value = lower + (upper - lower) / (1.0 + exp_neg_z)
  dvalue_dz = (upper - lower) * exp_neg_z / (1.0 + exp_neg_z) ** 2
  return value, dvalue_dz

def ss_from_f_obs(f_obs):
  """ (sin(theta)/lambda)^2 = d_star_sq/4, matching mmtbx.f_model.manager's
  internal self.ss exactly (f_model.py: self.ss = 1./flex.pow2(
  f_obs.d_spacings().data())/4.), recomputed here rather than reached for
  on the manager's internal .arrays.core.ss (no public accessor exists).
  """
  return 1.0 / flex.pow2(f_obs.d_spacings().data()) / 4.0

def initial_k_sol_b_sol(fmodel,
      k_sol_default=0.35, b_sol_default=46.0,
      k_sol_min=0.0, k_sol_max=0.6, b_sol_min=0.0, b_sol_max=150.0):
  """ Recover a starting (k_sol, b_sol) scalar pair from an fmodel's
  ALREADY-FITTED k_mask() array (e.g. right after the existing LS
  update_all_scales() call -- the design note sec. 8 call-0 bootstrap: an
  LS pre-fit places k_sol/b_sol in a reasonable neighbourhood before the
  LLGI-driven inner loop, sec. 7, takes over).

  Delegates directly to fmodel.k_sol_b_sol_from_k_mask() (mmtbx/f_model/
  f_model.py), bss's OWN summary of its own per-reflection k_mask() fit
  (a low-resolution-only Gaussian starting point, refined by a small
  local grid search on the sum-of-squares residual against k_sol*
  exp(-b_sol*ss), already clipped to [0, 0.6]/[0, 150] internally) --
  NOT reimplemented by hand here. An earlier version of this function
  did its own unweighted ordinary-least-squares fit of log(k_mask) =
  log(k_sol) - b_sol*ss directly against the raw per-reflection k_mask()
  array; on real data (2G38) this diverged badly (k_sol~2.1, b_sol~261,
  both then silently clipped to the k_sol_max/b_sol_max bounds) while
  fmodel.k_sol_b_sol_from_k_mask() on the SAME fmodel gave a sane
  (k_sol=0.43, b_sol=52.98) matching bss's own reported fit exactly --
  traced to bss's actual per-reflection k_mask() not being a clean
  k_sol*exp(-b_sol*ss) curve (see estimate_e_sigmaa_fixed_bulk_solvent's
  own docstring: "bss's bulk-solvent estimation is not guaranteed to
  already be an exact k_sol*exp(-b_sol*ss) form"), which an unweighted
  log-space OLS fit is highly sensitive to (near-zero mask values
  dominate the log-residual). That divergent starting point, once
  clipped to the phil bounds, was the actual root cause of the bulk-
  solvent LBFGS fit repeatedly "pinning at a phil boundary" documented
  elsewhere in this file: run_inner_loop's Stage 2 LBFGS was starting
  AT the boundary already (not degenerate at the optimum, and not
  reaching it via gradient ascent -- see the target-surface grid scan in
  doc/llgi_target_design.md's bulk-solvent diagnostic note, which found
  a well-defined interior minimum near (0.35-0.40, 46-80) that was
  simply never reached because LBFGS started outside it and the bounded-
  sigmoid reparameterisation's gradient vanishes near the tails, per
  _bounded's own docstring). Fixed by reusing the existing, correct
  accessor instead of re-deriving a fragile approximation to it.

  Returns (k_sol, b_sol), both plain floats, already clipped to
  [0, 0.6]/[0, 150] by k_sol_b_sol_from_k_mask() itself; the k_sol_min/
  max/b_sol_min/max arguments are kept for backward API compatibility
  (callers already pass params.k_sol_min etc.) and applied as an
  additional clip in case a caller's phil bounds are narrower than
  k_sol_b_sol_from_k_mask()'s own fixed [0, 0.6]/[0, 150] range.
  """
  k_masks = fmodel.k_masks()
  assert len(k_masks) == 1, (
    "initial_k_sol_b_sol: only a single bulk-solvent mask shell is "
    "supported (design note sec. 9); got %d." % len(k_masks))
  k_sol, b_sol = fmodel.k_sol_b_sol_from_k_mask()
  if(k_sol is None or b_sol is None):
    return k_sol_default, b_sol_default
  k_sol = min(max(float(k_sol), k_sol_min), k_sol_max)
  b_sol = min(max(float(b_sol), b_sol_min), b_sol_max)
  return k_sol, b_sol

def k_mask_and_gradients(ss, k_sol, b_sol):
  """ k_mask(ss) = k_sol*exp(-b_sol*ss) (ext.k_mask's own formula,
  mmtbx/f_model/f_model.h's f_b_exp_one_h/k_mask), plus its two partial
  derivatives, needed by the Stage-2 (bulk-solvent) chain rule (design
  note sec. 6):
    d(k_mask)/d(k_sol) = exp(-b_sol*ss) = k_mask/k_sol  (k_sol != 0)
    d(k_mask)/d(b_sol) = -ss*k_sol*exp(-b_sol*ss) = -ss*k_mask

  Returns a group_args with .k_mask, .d_by_dk_sol, .d_by_db_sol (all
  flex.double, same order as ss). Reuses ext.k_mask itself (rather than
  reimplementing exp(-b_sol*ss) by hand) so the value exactly matches
  what the rest of mmtbx computes, including any future changes to that
  function's numerical details (e.g. exponent clamping).
  """
  k_mask = f_model_ext.k_mask(ss, k_sol, b_sol)
  d_by_db_sol = -ss * k_mask
  if(k_sol != 0):
    d_by_dk_sol = k_mask / k_sol
  else:
    d_by_dk_sol = flex.exp(-b_sol * ss)
  return group_args(
    k_mask=k_mask, d_by_dk_sol=d_by_dk_sol, d_by_db_sol=d_by_db_sol)

def bulk_solvent_target_and_gradients(
      e_eff, selection, dobs, sigmaa, centric_flags,
      f_calc, f_mask, f_part1, f_part2, k_isotropic,
      epsilons, d_star_sq, ss, k_sol, b_sol,
      n_sigmap_nodes=15, auto_kernel_number=50):
  """ E-scale LLGI target and its gradient w.r.t. (k_sol, b_sol), for the
  Stage-2 (bulk-solvent) fit (design note sec. 6), with sigmaA(d) fixed
  (already evaluated per reflection -- the whole point of this stage is
  to hold it constant, sec. 7).

  Single mask shell only (design note sec. 9). f_calc, f_mask, f_part1,
  f_part2, k_isotropic are the already-scaled-except-bulk-solvent
  components (fmodel.f_calc(), fmodel.f_masks()[0], fmodel.f_part1(),
  fmodel.f_part2(), fmodel.k_isotropic() respectively) -- held fixed
  within one gradient evaluation, exactly like SigmaP (see below).

  Chain rule (design note sec. 6):
    d(target)/d(k_sol) = sum_i [ d_target_over_demodel_i *
      Re( conj(Emodel_i)/|Emodel_i| * d(Emodel_i)/d(k_sol) ) ]
  and likewise for b_sol, where:
    d(Emodel_i)/d(param) = (1/sqrt(EPS_i*SigmaP)) *
      d(f_model_no_aniso_scale_i)/d(param)
    d(fmnas_i)/d(param) = k_isotropic_i * f_mask_i * d(k_mask_i)/d(param)
  SigmaP is treated as locally fixed WITHIN this one gradient evaluation
  (not differentiated through -- differentiating the kernel-weighted
  average and Chebyshev fit analytically would be a substantially bigger
  undertaking), but this function itself recomputes SigmaP fresh from
  the k_sol/b_sol passed in every single time it is called -- so a
  caller (e.g. bulk_solvent_target_evaluator below) that calls this once
  per LBFGS step always sees an up-to-date target value at its current
  point, even though the gradient direction is the cheaper fixed-SigmaP
  approximation. This distinction matters: verified by finite difference
  that the fixed-SigmaP gradient (what this function returns) and the
  TRUE total derivative (letting SigmaP respond to k_sol/b_sol too) can
  differ by roughly an order of magnitude on real data -- SigmaP is a
  smoothed function of |f_model_no_aniso_scale|^2, which obviously
  depends on k_sol/b_sol directly, so this is a real approximation, not
  a negligible numerical nuance. The fixed-SigmaP finite-difference
  check (tst_llgi_e_bulk_solvent.py's exercise_bulk_solvent_gradient_
  matches_fixed_sigma_p_finite_difference) verifies THIS function's
  gradient is internally consistent with THIS function's own target
  (agreement to ~1e-6 relative), not with the true total derivative.

  Returns a group_args with .target (float) and .gradients (numpy array,
  [d_target_by_dk_sol, d_target_by_db_sol]).
  """
  import numpy as np
  k_mask_result = k_mask_and_gradients(ss, k_sol, b_sol)
  bulk = k_mask_result.k_mask * f_mask
  fmnas = k_isotropic * (f_calc + bulk + f_part1 + f_part2)
  e_model_result = build_e_model(
    fmnas, epsilons, d_star_sq,
    n_sigmap_nodes=n_sigmap_nodes, auto_kernel_number=auto_kernel_number)
  e_model = e_model_result.e_model
  sigma_p = e_model_result.sigma_p
  e_model_abs = flex.abs(e_model)

  result = xray_ext.llgi_e_emodel_target_and_gradients(
    e_eff=e_eff,
    selection=selection,
    e_model=e_model_abs,
    dobs=dobs,
    sigmaa=sigmaa,
    centric_flags=centric_flags)
  d_target_by_demodel = np.array(result.d_target_by_demodel())

  # d(fmnas)/d(param) = k_isotropic * f_mask * d(k_mask)/d(param)
  d_fmnas_by_dk_sol = k_isotropic * f_mask * k_mask_result.d_by_dk_sol
  d_fmnas_by_db_sol = k_isotropic * f_mask * k_mask_result.d_by_db_sol
  inv_denom = 1.0 / flex.sqrt(epsilons * sigma_p)
  d_emodel_by_dk_sol = d_fmnas_by_dk_sol * inv_denom
  d_emodel_by_db_sol = d_fmnas_by_db_sol * inv_denom

  # Real-to-complex phase projection, matching llgi.h's
  # d_target_one_h_over_fc pattern: Re( conj(E)/|E| * dE/dparam ). flex
  # does not support complex_double / double directly, hence the
  # reciprocal-and-multiply below (as elsewhere in this module).
  e_model_abs_safe = flex.double([
    v if v > 0 else 1.0 for v in e_model_abs])
  phase = flex.conj(e_model) * (1.0 / e_model_abs_safe)
  proj_k_sol = flex.real(phase * d_emodel_by_dk_sol)
  proj_b_sol = flex.real(phase * d_emodel_by_db_sol)

  # Only selected reflections (and only where Emodel > 0, matching the
  # C++ target's own dobs/sigmaa/eeff/emodel > 0 guards, which silently
  # zero d_target_by_demodel there already) contribute.
  d_target_by_dk_sol = float(np.dot(d_target_by_demodel,
    proj_k_sol.as_numpy_array()))
  d_target_by_db_sol = float(np.dot(d_target_by_demodel,
    proj_b_sol.as_numpy_array()))
  return group_args(
    target=result.target(),
    gradients=np.array([d_target_by_dk_sol, d_target_by_db_sol]))

class bulk_solvent_target_evaluator(object):
  """ scitbx.lbfgs target evaluator optimising (k_sol, b_sol) against the
  E-scale LLGI target, summed over the working set (all reflections
  excluding R-free -- design note sec. 6), with sigmaA(d) held fixed.
  x = [z_k_sol, z_b_sol], unconstrained; actual (k_sol, b_sol) obtained
  via _bounded (design note doesn't specify a reparameterisation for
  bulk solvent explicitly, but the same bounded-sigmoid convention used
  for sigmaA is reused here for consistency and to keep LBFGS well-
  behaved against mmtbx's usual k_sol/b_sol ranges).
  """

  def __init__(self,
        e_eff, working_selection, dobs, sigmaa, centric_flags,
        f_calc, f_mask, f_part1, f_part2, k_isotropic,
        epsilons, d_star_sq, ss,
        k_sol_start, b_sol_start,
        k_sol_min=0.0, k_sol_max=0.6, b_sol_min=0.0, b_sol_max=150.0,
        n_sigmap_nodes=15, auto_kernel_number=50,
        max_iterations=50):
    self.e_eff = e_eff
    self.working_selection = working_selection
    self.dobs = dobs
    self.sigmaa = sigmaa
    self.centric_flags = centric_flags
    self.f_calc = f_calc
    self.f_mask = f_mask
    self.f_part1 = f_part1
    self.f_part2 = f_part2
    self.k_isotropic = k_isotropic
    self.epsilons = epsilons
    self.d_star_sq = d_star_sq
    self.ss = ss
    self.k_sol_min = k_sol_min
    self.k_sol_max = k_sol_max
    self.b_sol_min = b_sol_min
    self.b_sol_max = b_sol_max
    self.n_sigmap_nodes = n_sigmap_nodes
    self.auto_kernel_number = auto_kernel_number
    import numpy as np
    def _inv_bounded(value, lower, upper):
      # Inverse of _bounded: z such that _bounded(z) == value.
      frac = (value - lower) / (upper - lower)
      frac = min(max(frac, 1.e-6), 1.0 - 1.e-6)
      return float(np.log(frac / (1.0 - frac)))
    self.x = flex.double([
      _inv_bounded(k_sol_start, k_sol_min, k_sol_max),
      _inv_bounded(b_sol_start, b_sol_min, b_sol_max)])
    self.final_target = None
    term_parameters = scitbx.lbfgs.termination_parameters(
      max_iterations=max_iterations)
    exception_handling_parameters = scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound=True,
      ignore_line_search_failed_step_at_upper_bound=True)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=term_parameters,
      exception_handling_params=exception_handling_parameters)

  def _current_k_sol_b_sol(self):
    import numpy as np
    x_np = np.array(self.x)
    k_sol, dk_sol_dz = _bounded(x_np[0], self.k_sol_min, self.k_sol_max)
    b_sol, db_sol_dz = _bounded(x_np[1], self.b_sol_min, self.b_sol_max)
    return (float(k_sol), float(dk_sol_dz)), (float(b_sol), float(db_sol_dz))

  def compute_functional_and_gradients(self):
    (k_sol, dk_sol_dz), (b_sol, db_sol_dz) = self._current_k_sol_b_sol()
    result = bulk_solvent_target_and_gradients(
      e_eff=self.e_eff,
      selection=self.working_selection,
      dobs=self.dobs,
      sigmaa=self.sigmaa,
      centric_flags=self.centric_flags,
      f_calc=self.f_calc,
      f_mask=self.f_mask,
      f_part1=self.f_part1,
      f_part2=self.f_part2,
      k_isotropic=self.k_isotropic,
      epsilons=self.epsilons,
      d_star_sq=self.d_star_sq,
      ss=self.ss,
      k_sol=k_sol,
      b_sol=b_sol,
      n_sigmap_nodes=self.n_sigmap_nodes,
      auto_kernel_number=self.auto_kernel_number)
    self.final_target = result.target
    g_k_sol = result.gradients[0] * dk_sol_dz
    g_b_sol = result.gradients[1] * db_sol_dz
    return result.target, flex.double([g_k_sol, g_b_sol])

  def k_sol_b_sol(self):
    (k_sol, _), (b_sol, _) = self._current_k_sol_b_sol()
    return k_sol, b_sol

class e_sigmaa_target_evaluator(object):
  """ scitbx.lbfgs target evaluator optimising the B-spline coefficients
  of the E-scale sigmaA(resolution) curve against the E-scale LLGI
  target, summed over the R-free/test set only (design note sec. 5),
  with Emodel (hence the bulk-solvent model) held fixed. Directly
  mirrors mmtbx.refinement.llgi_sigmaa.llgi_sigmaa_target_evaluator
  (same B-spline-over-sigmoid parameterisation, same LBFGS setup); the
  only change is which target functor supplies the LLG value/gradient
  (ext.llgi_e_sigmaa_target_and_gradients here, vs
  ext.llgi_sigmaa_scatfrac_target_and_gradients there).
  """

  def __init__(self,
        e_eff, test_selection, e_model, dobs, centric_flags,
        sigmaa_design, n_sigmaa_coeffs, max_iterations=100,
        curvature_weight=0.0, spline_degree=3):
    self.e_eff = e_eff
    self.test_selection = test_selection
    self.e_model = e_model
    self.dobs = dobs
    self.centric_flags = centric_flags
    self.sigmaa_design = sigmaa_design  # numpy array, (n_refl, n_coeffs)
    self.n_sigmaa_coeffs = n_sigmaa_coeffs
    self.curvature_weight = curvature_weight
    # Only needed by evaluate_at() (re-evaluating the fitted curve at
    # NEW d_star_sq values, e.g. missing reflections for map-coefficient
    # fill-missing) -- NOT used anywhere in the fit itself, which only
    # ever consults the already-built sigmaa_design above; kept as a
    # plain default-valued constructor argument (not re-derived from
    # sigmaa_design, which has no record of what degree built it) so
    # existing callers that don't pass it keep working unchanged.
    self._spline_degree = spline_degree
    # Unconstrained starting point: z=0 maps (via the sigmoid) to
    # sigmaA=0.5, a neutral starting guess (matches llgi_sigmaa's).
    self.x = flex.double(n_sigmaa_coeffs, 0.0)
    term_parameters = scitbx.lbfgs.termination_parameters(
      max_iterations=max_iterations)
    exception_handling_parameters = scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound=True,
      ignore_line_search_failed_step_at_upper_bound=True)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=term_parameters,
      exception_handling_params=exception_handling_parameters)

  def _current_sigmaa(self):
    import numpy as np
    coeffs = np.array(self.x)
    z = self.sigmaa_design.dot(coeffs)
    sigmaa, dsigmaa_dz = _sigmoid(z)
    return sigmaa, dsigmaa_dz

  def compute_functional_and_gradients(self):
    import numpy as np
    sigmaa, dsigmaa_dz = self._current_sigmaa()
    result = xray_ext.llgi_e_sigmaa_target_and_gradients(
      e_eff=self.e_eff,
      selection=self.test_selection,
      e_model=self.e_model,
      dobs=self.dobs,
      sigmaa=flex.double(sigmaa),
      centric_flags=self.centric_flags)
    f = result.target()
    d_target_by_dsigmaa = np.array(result.d_target_by_dsigmaa())
    g = self.sigmaa_design.T.dot(d_target_by_dsigmaa * dsigmaa_dz)
    penalty, penalty_grad = _spline_curvature_penalty_and_gradient(
      np.array(self.x), self.curvature_weight)
    f += penalty
    g = g + penalty_grad
    return f, flex.double(g)

  def sigmaa(self):
    sigmaa, _ = self._current_sigmaa()
    return flex.double(sigmaa)

  def evaluate_at(self, d_star_sq, x_range):
    """ Evaluate the ALREADY-FITTED sigmaA(d) curve (this evaluator's own
    converged B-spline coefficients, self.x) at an arbitrary new set of
    d_star_sq values -- e.g. a missing-reflection index set, for map-
    coefficient fill-missing (mmtbx.map_tools.model_missing_reflections_
    llgi) -- rather than the design matrix this evaluator was actually
    fit against (self.sigmaa_design, built from the OBSERVED reflection
    set's own d_star_sq).

    x_range MUST be the same (x_min, x_max) pair the fit's own design
    matrix was built with (see _b_spline_design_matrix's own docstring
    on why this must be passed explicitly and consistently) -- the
    caller (estimate_e_sigmaa) records this as .x_range in its own
    result for exactly this purpose.

    Points outside x_range are CLAMPED to the nearest endpoint's fitted
    value, not extrapolated: _b_spline_design_matrix builds its design
    matrix with scipy's extrapolate=False, so a genuinely out-of-range
    point would otherwise silently get an all-zero design row (and
    hence sigmaA=_sigmoid(0)=0.5, a meaningless default, not a curve
    value) -- clamping avoids that failure mode for e.g. a missing
    reflection at lower resolution than anything in the observed set
    (a plausible case for systematic absences / detector gaps).

    Returns a flex.double, one sigmaA value per input d_star_sq.
    """
    import numpy as np
    d_star_sq_np = np.asarray(d_star_sq, dtype=float)
    x_min, x_max = x_range
    d_star_sq_clamped = np.clip(d_star_sq_np, x_min, x_max)
    design = _b_spline_design_matrix(
      d_star_sq_clamped, self.n_sigmaa_coeffs, self._spline_degree,
      x_range=x_range)
    coeffs = np.array(self.x)
    z = design.dot(coeffs)
    sigmaa, _ = _sigmoid(z)
    return flex.double(sigmaa)

def estimate_e_sigmaa(e_eff, r_free_flags, e_model, dobs, centric_flags,
      d_star_sq, n_coeffs=8, spline_degree=3, max_iterations=100,
      curvature_weight=0.0):
  """ Fit the E-scale sigmaA(resolution) curve against the E-scale LLGI
  target, restricted to the R-free/test set (design note sec. 5), with
  Emodel (i.e. the current bulk-solvent model) held fixed. Evaluates the
  fitted curve at every reflection.

  Returns a group_args with .sigmaa (flex.double, one value per input
  reflection), .target (final fitted LLGI target value on the test set,
  for diagnostics/logging), .x_range (the (d_star_sq_min, d_star_sq_max)
  the fit's B-spline design matrix was actually built against), and
  .evaluate_at (a bound method, evaluator.evaluate_at, for re-evaluating
  this SAME fitted curve at NEW d_star_sq values -- e.g. mmtbx.map_tools
  .model_missing_reflections_llgi's own missing-reflection fill; pass
  x_range=result.x_range to it, exactly as documented on evaluate_at's
  own docstring).
  """
  n_refl = e_eff.size()
  assert r_free_flags.size() == n_refl
  assert e_model.size() == n_refl
  assert dobs.size() == n_refl
  assert centric_flags.size() == n_refl
  assert d_star_sq.size() == n_refl
  n_test = r_free_flags.count(True)
  if(n_test == 0):
    raise RuntimeError(
      "No R-free/test-set reflections available for the E-scale LLGI "
      "sigmaA fit.")
  d_star_sq_np = d_star_sq.as_numpy_array()
  x_range = (float(d_star_sq_np.min()), float(d_star_sq_np.max()))
  sigmaa_design = _b_spline_design_matrix(
    d_star_sq_np, n_coeffs, spline_degree, x_range=x_range)
  evaluator = e_sigmaa_target_evaluator(
    e_eff=e_eff,
    test_selection=r_free_flags,
    e_model=e_model,
    dobs=dobs,
    centric_flags=centric_flags,
    sigmaa_design=sigmaa_design,
    n_sigmaa_coeffs=n_coeffs,
    max_iterations=max_iterations,
    curvature_weight=curvature_weight,
    spline_degree=spline_degree)
  sigmaa = evaluator.sigmaa()
  final_result = xray_ext.llgi_e_sigmaa_target_and_gradients(
    e_eff=e_eff, selection=r_free_flags, e_model=e_model, dobs=dobs,
    sigmaa=sigmaa, centric_flags=centric_flags)
  return group_args(
    sigmaa=sigmaa, target=final_result.target(), x_range=x_range,
    evaluate_at=evaluator.evaluate_at)

def estimate_e_sigmaa_fixed_bulk_solvent(
      fmodel, dobs, feff, resn, params=None):
  """ Hybrid mode (params.fix_bulk_solvent_from_ls): fit the E-scale
  sigmaA(resolution) curve (Stage 1 only) against fmodel's bulk-solvent
  model EXACTLY as it currently stands -- i.e. bss's own least-squares
  (k_sol, b_sol) (or, more precisely, its own per-reflection k_mask
  array, whatever shape that took: bss's bulk-solvent estimation is not
  guaranteed to already be an exact k_sol*exp(-b_sol*ss) form, so this
  reads fmodel.k_masks()[0] directly via f_model_no_aniso_scale rather
  than going through a (k_sol, b_sol) scalar re-fit the way
  initial_k_sol_b_sol/run_inner_loop's Stage 2 does -- there is no
  approximation of bss's result here, unlike the normal inner loop's
  bootstrap). Stage 2 (the LLGI bulk-solvent fit) never runs at all, and
  fmodel's k_mask is never modified -- see the design note sec. 6's "LS-
  fixed bulk solvent, LLGI-only sigmaA" note for why: repeated testing
  (mmtbx.regression.tst_llgi_e_bulk_solvent and real-data macrocycle
  runs) found Stage 2 driving (k_sol, b_sol) to one of two competing
  solutions, one or both pinned at a phil boundary, with a fitted-curve-
  smoothness restraint weight that helps at the first macrocycle but
  stops mattering by the second -- this mode sidesteps that failure mode
  entirely by never letting LLGI touch bulk solvent.

  fmodel: an mmtbx.f_model.manager, already scaled (matching
  run_inner_loop's own precondition -- e.g. straight after
  update_all_scales()/bss).

  dobs, feff, resn: as run_inner_loop (already matched to fmodel.f_obs()'s
  current index set -- see run_inner_loop's docstring for the same
  index-set-drift hazard, unchanged here).

  params: extracted llgi_e_bulk_solvent_params phil, or None for
  defaults. k_sol_min/max, b_sol_min/max, bulk_solvent_max_iterations,
  max_inner_iterations, convergence_tolerance are all unused here (no
  Stage 2, no iteration).

  Returns a group_args with the same shape as run_inner_loop's result
  (.sigmaa, .k_sol, .b_sol, .n_iterations, .converged, .history) so
  callers (mmtbx.f_model.manager.update_llgi_e_bulk_solvent) do not need
  to special-case this mode -- .k_sol/.b_sol are a log-linear POINT
  ESTIMATE of bss's actual k_mask (via initial_k_sol_b_sol, for logging/
  diagnostics only -- NOT what was actually used to build Emodel, which
  is bss's raw per-reflection k_mask), .n_iterations=1, .converged=True,
  .history has a single entry.
  """
  if(params is None):
    params = llgi_e_bulk_solvent_params.extract()
  f_obs = fmodel.f_obs()
  r_free_flags = fmodel.r_free_flags().data()
  centric_flags = f_obs.centric_flags().data()
  epsilons = f_obs.epsilons().data().as_double()
  d_star_sq = f_obs.d_star_sq().data()

  e_eff = build_e_eff(feff, resn)
  fmnas = f_model_no_aniso_scale(fmodel).data()
  e_model = build_e_model(
    fmnas, epsilons, d_star_sq,
    n_sigmap_nodes=params.n_sigmap_nodes,
    auto_kernel_number=params.auto_kernel_number).e_model

  sigmaa_result = estimate_e_sigmaa(
    e_eff=e_eff, r_free_flags=r_free_flags,
    e_model=flex.abs(e_model), dobs=dobs,
    centric_flags=centric_flags, d_star_sq=d_star_sq,
    n_coeffs=params.n_sigmaa_coeffs, spline_degree=params.spline_degree,
    max_iterations=params.sigmaa_max_iterations,
    curvature_weight=params.sigmaa_curvature_weight)

  # Point estimate only, for logging -- see docstring above; not used to
  # build e_model (bss's raw k_mask was used directly instead).
  k_sol, b_sol = initial_k_sol_b_sol(
    fmodel,
    k_sol_min=params.k_sol_min, k_sol_max=params.k_sol_max,
    b_sol_min=params.b_sol_min, b_sol_max=params.b_sol_max)

  return group_args(
    sigmaa=sigmaa_result.sigmaa,
    k_sol=k_sol, b_sol=b_sol,
    n_iterations=1,
    converged=True,
    history=[group_args(
      sigmaa_target=sigmaa_result.target, bs_target=None,
      k_sol=k_sol, b_sol=b_sol)])

def estimate_sigmaa_e_then_scatfrac_f(
      fmodel, dobs, feff, resn, e_params=None, scatfrac_params=None):
  """ Break the F-scale sigmaA/ScatFrac non-identifiability by fitting
  sigmaA FIRST against the E-scale LLGI target, then ScatFrac against the
  F-scale LLGI target with sigmaA fixed -- the alternative to the
  historical scheme (ScatFrac by moment estimator, then sigmaA by F-scale
  LLGI), selected by mmtbx.refinement.llgi_sigmaa.
  llgi_sigmaa_scatfrac_params.estimate_scatfrac_by_likelihood.

  Why this ordering avoids the degeneracy without needing a moment
  estimator: the F-scale target only sees sigmaA and ScatFrac through the
  combination D = Dobs*sigmaA/sqrt(ScatFrac), so a joint fit of both
  against it is degenerate. The E-scale target (llgi_e.h/estimate_e_
  sigmaa) has NO ScatFrac term at all -- Emodel is normalised by
  sqrt(EPS*SigmaP), not by ScatFrac -- so fitting sigmaA there first is
  not a degenerate problem, and the F-scale ScatFrac fit that follows has
  only one free curve left.

  Step 1 (sigmaA, E-scale, R-free only): delegates to estimate_e_sigmaa_
  fixed_bulk_solvent, i.e. bulk solvent is whatever fmodel's CURRENT
  state already is (bss's LS fit, per the "LS-fixed bulk solvent" default
  -- see that function's docstring) -- this function does not touch bulk
  solvent itself, consistent with keeping this a sigmaA/ScatFrac
  estimator only, not a bulk-solvent one.

  Step 2 (ScatFrac, F-scale, working set = all minus R-free): delegates
  to mmtbx.refinement.llgi_sigmaa.estimate_llgi_scatfrac_likelihood,
  fixing sigmaA at Step 1's result. Uses fmodel.f_model() (bulk-solvent-
  and scale-corrected), matching update_llgi_sigmaa_scatfrac's own
  convention and rationale (raw f_calc() lacks the k_isotropic correction
  Feff/Resn/f_model() already reflect -- see that method's docstring for
  the real-data bug this was traced to previously). ScatFrac is NOT
  bounded above by 1: Feff need not be on absolute scale (its scaling
  assumes 50% solvent content by default), so ScatFrac can genuinely
  exceed 1 -- see llgi_scatfrac_target_evaluator's docstring.

  fmodel: an mmtbx.f_model.manager, already scaled (e.g. straight after
  update_all_scales()/bss).

  dobs, feff, resn: nacelle DOBS/FEFF/RESN columns, already matched to
  fmodel.f_obs()'s CURRENT index set (same index-set-drift hazard as
  run_inner_loop/estimate_e_sigmaa_fixed_bulk_solvent -- see those
  functions' docstrings).

  e_params: extracted llgi_e_bulk_solvent_params phil, or None for
  defaults (passed straight through to estimate_e_sigmaa_fixed_bulk_
  solvent for Step 1).

  scatfrac_params: extracted llgi_sigmaa_scatfrac_params phil, or None
  for defaults (passed straight through to estimate_llgi_scatfrac_
  likelihood for Step 2; only n_scatfrac_coeffs, spline_degree,
  n_scatfrac_bins, max_iterations are used there).

  Returns a group_args with .sigmaa, .scatfrac (flex.double, one value
  per input reflection, evaluated at every reflection), .target (the
  Step 2 ScatFrac fit's final LLGI target value on the working set, for
  diagnostics/logging -- NOT directly comparable to update_llgi_sigmaa_
  scatfrac's historical .target, which is evaluated on R-free), and
  .scatfrac_inf/.b_scatfrac (floats, only set when scatfrac_params.
  scatfrac_model == "b_factor"; None otherwise -- forwarded from
  estimate_llgi_scatfrac_likelihood's own result for logging/
  diagnostics, see llgi_scatfrac_b_factor_target_evaluator's docstring
  for why these two values should not be over-interpreted on their own).
  """
  import mmtbx.refinement.llgi_sigmaa as llgi_sigmaa
  if(e_params is None):
    e_params = llgi_e_bulk_solvent_params.extract()
  if(scatfrac_params is None):
    scatfrac_params = llgi_sigmaa.llgi_sigmaa_scatfrac_params.extract()

  sigmaa_result = estimate_e_sigmaa_fixed_bulk_solvent(
    fmodel, dobs=dobs, feff=feff, resn=resn, params=e_params)
  sigmaa = sigmaa_result.sigmaa

  f_obs = fmodel.f_obs()
  r_free_flags = fmodel.r_free_flags().data()
  working_selection = ~r_free_flags
  llgi_data = fmodel.llgi_data()
  teps = llgi_data.teps.data()
  scatfrac_result = llgi_sigmaa.estimate_llgi_scatfrac_likelihood(
    f_eff=feff,
    working_selection=working_selection,
    f_calc=fmodel.f_model().data(),
    dobs=dobs,
    sigmaa=sigmaa,
    teps=teps,
    resn=resn,
    centric_flags=f_obs.centric_flags().data(),
    d_star_sq=f_obs.d_star_sq().data(),
    scale_factor=fmodel.scale_ml_wrapper(),
    params=scatfrac_params)

  return group_args(
    sigmaa=sigmaa, scatfrac=scatfrac_result.scatfrac,
    target=scatfrac_result.target,
    scatfrac_inf=getattr(scatfrac_result, "scatfrac_inf", None),
    b_scatfrac=getattr(scatfrac_result, "b_scatfrac", None))

def run_inner_loop(fmodel, dobs, feff, resn, params=None, log=None):
  """ The full Stage-1/Stage-2 inner loop (design note sec. 7): alternate
  fitting sigmaA(d) (E-scale, R-free only, sec. 5) and the bulk-solvent
  model (k_sol, b_sol; E-scale, working set only, sec. 6) to local
  convergence within one macrocycle, holding the other fixed at each
  step.

  fmodel: an mmtbx.f_model.manager, already scaled (k_isotropic/
  k_anisotropic and an initial bulk-solvent estimate should already be
  in place -- e.g. from the existing LS update_all_scales(), matching
  the call-0 bootstrap in the design note sec. 8. Only overall scale/
  anisotropy are left untouched by this function (sec. 9): only k_sol/
  b_sol (via fmodel.update(k_mask=...)) are modified here.

  dobs, feff, resn: nacelle DOBS/FEFF/RESN columns, already matched to
  fmodel.f_obs()'s CURRENT index set (see phenix.refinement.llgi_data.
  get_llgi_data) -- specifically, the index set AFTER any outlier
  removal fmodel's own update_all_scales() may have already performed
  (update_all_scales(remove_outliers=True), its default, can shrink
  f_obs -- observed directly while testing this function: a 16480-
  reflection f_obs became 16477 after one update_all_scales() call).
  Passing dobs/feff/resn sized against a PRE-outlier-removal f_obs
  raises AssertionError here rather than silently misaligning HKLs; the
  same index-set-drift hazard the earlier llgi_data/select() fix
  addressed for the F-scale sigmaa/scatfrac estimator (see doc/
  llgi_target_design.md sec. 5's "llgi_data dropped by manager.select()"
  bug and mmtbx.f_model.manager.set_llgi_data/llgi_data). Eeff = Feff/
  RESN is built once, up front (sec. 4: the experimental data does not
  change within this loop).

  params: extracted llgi_e_bulk_solvent_params phil, or None for
  defaults.

  Returns a group_args with:
    .sigmaa       -- final fitted sigmaA(d) (flex.double)
    .k_sol, .b_sol -- final bulk-solvent parameters (floats)
    .n_iterations -- number of Stage-1/Stage-2 alternations actually run
    .converged    -- bool
    .history      -- list of group_args(sigmaa_target=..., bs_target=...,
                      k_sol=..., b_sol=...), one entry per iteration, for
                      diagnostics/logging.
  Mutates fmodel in place (updates its k_mask via fmodel.update()) to
  reflect the final bulk-solvent fit; sigmaA itself is not stored on
  fmodel (it is not one of fmodel's own scale parameters), just returned.

  If params.fix_bulk_solvent_from_ls is True, this delegates entirely to
  estimate_e_sigmaa_fixed_bulk_solvent instead (Stage 2/bulk-solvent
  fitting skipped, fmodel's k_mask left untouched -- see that function's
  docstring); the rest of this docstring describes the default
  (Stage-1/Stage-2) behaviour only.
  """
  if(params is None):
    params = llgi_e_bulk_solvent_params.extract()
  if(params.fix_bulk_solvent_from_ls):
    return estimate_e_sigmaa_fixed_bulk_solvent(
      fmodel, dobs=dobs, feff=feff, resn=resn, params=params)
  f_obs = fmodel.f_obs()
  r_free_flags = fmodel.r_free_flags().data()
  centric_flags = f_obs.centric_flags().data()
  epsilons = f_obs.epsilons().data().as_double()
  d_star_sq = f_obs.d_star_sq().data()
  ss = ss_from_f_obs(f_obs)
  working_selection = ~r_free_flags

  e_eff = build_e_eff(feff, resn)

  k_sol, b_sol = initial_k_sol_b_sol(
    fmodel,
    k_sol_min=params.k_sol_min, k_sol_max=params.k_sol_max,
    b_sol_min=params.b_sol_min, b_sol_max=params.b_sol_max)

  history = []
  prev_bs_target = None
  converged = False
  n_iterations = 0

  for iteration in range(params.max_inner_iterations):
    n_iterations = iteration + 1
    f_calc = fmodel.f_calc().data()
    f_masks = fmodel.f_masks()
    assert len(f_masks) == 1, (
      "run_inner_loop: only a single bulk-solvent mask shell is "
      "supported (design note sec. 9); got %d." % len(f_masks))
    f_mask = f_masks[0].data()
    f_part1 = fmodel.f_part1().data()
    f_part2 = fmodel.f_part2().data()
    k_isotropic = fmodel.k_isotropic()

    # Stage 1: fit sigmaA(d) against LLG on R-free, Emodel (hence bulk
    # solvent) held fixed at this iteration's current k_sol/b_sol.
    k_mask_current = k_mask_and_gradients(ss, k_sol, b_sol).k_mask
    fmnas_current = k_isotropic * (
      f_calc + k_mask_current * f_mask + f_part1 + f_part2)
    e_model_current = build_e_model(
      fmnas_current, epsilons, d_star_sq,
      n_sigmap_nodes=params.n_sigmap_nodes,
      auto_kernel_number=params.auto_kernel_number).e_model
    sigmaa_result = estimate_e_sigmaa(
      e_eff=e_eff, r_free_flags=r_free_flags,
      e_model=flex.abs(e_model_current), dobs=dobs,
      centric_flags=centric_flags, d_star_sq=d_star_sq,
      n_coeffs=params.n_sigmaa_coeffs, spline_degree=params.spline_degree,
      max_iterations=params.sigmaa_max_iterations,
      curvature_weight=params.sigmaa_curvature_weight)
    sigmaa = sigmaa_result.sigmaa

    # Stage 2: fit (k_sol, b_sol) against LLG on the working set,
    # sigmaA(d) held fixed at the value just obtained above.
    bs_evaluator = bulk_solvent_target_evaluator(
      e_eff=e_eff, working_selection=working_selection, dobs=dobs,
      sigmaa=sigmaa, centric_flags=centric_flags,
      f_calc=f_calc, f_mask=f_mask, f_part1=f_part1, f_part2=f_part2,
      k_isotropic=k_isotropic, epsilons=epsilons, d_star_sq=d_star_sq,
      ss=ss, k_sol_start=k_sol, b_sol_start=b_sol,
      k_sol_min=params.k_sol_min, k_sol_max=params.k_sol_max,
      b_sol_min=params.b_sol_min, b_sol_max=params.b_sol_max,
      n_sigmap_nodes=params.n_sigmap_nodes,
      auto_kernel_number=params.auto_kernel_number,
      max_iterations=params.bulk_solvent_max_iterations)
    k_sol, b_sol = bs_evaluator.k_sol_b_sol()
    bs_target = bs_evaluator.final_target

    if(log is not None):
      print(
        "LLGI E-scale inner loop iter %d: sigmaA target=%.6f  "
        "bulk-solvent target=%.6f  k_sol=%.4f  b_sol=%.2f" % (
          n_iterations, sigmaa_result.target, bs_target, k_sol, b_sol),
        file=log)
    history.append(group_args(
      sigmaa_target=sigmaa_result.target, bs_target=bs_target,
      k_sol=k_sol, b_sol=b_sol))

    if(prev_bs_target is not None and
       abs(bs_target - prev_bs_target) < params.convergence_tolerance):
      converged = True
      break
    prev_bs_target = bs_target

  # Push the final (k_sol, b_sol) back onto fmodel via its k_mask, per
  # the design note's chosen mechanism (fmodel.update(k_mask=...) --
  # mmtbx.f_model.manager has no direct k_sol/b_sol setter; see
  # k_mask_and_gradients/initial_k_sol_b_sol's docstrings).
  final_k_mask = k_mask_and_gradients(ss, k_sol, b_sol).k_mask
  fmodel.update(k_mask=[final_k_mask])

  # Final sigmaA fit at the converged bulk-solvent model, so the
  # returned .sigmaa reflects the LAST bulk-solvent state exactly (the
  # loop above may have exited right after a Stage-2 update, whose
  # matching Stage-1 refit would be the next iteration's first half).
  f_calc = fmodel.f_calc().data()
  f_mask = fmodel.f_masks()[0].data()
  f_part1 = fmodel.f_part1().data()
  f_part2 = fmodel.f_part2().data()
  k_isotropic = fmodel.k_isotropic()
  fmnas_final = k_isotropic * (f_calc + final_k_mask * f_mask
    + f_part1 + f_part2)
  e_model_final = build_e_model(
    fmnas_final, epsilons, d_star_sq,
    n_sigmap_nodes=params.n_sigmap_nodes,
    auto_kernel_number=params.auto_kernel_number).e_model
  final_sigmaa_result = estimate_e_sigmaa(
    e_eff=e_eff, r_free_flags=r_free_flags,
    e_model=flex.abs(e_model_final), dobs=dobs,
    centric_flags=centric_flags, d_star_sq=d_star_sq,
    n_coeffs=params.n_sigmaa_coeffs, spline_degree=params.spline_degree,
    max_iterations=params.sigmaa_max_iterations,
    curvature_weight=params.sigmaa_curvature_weight)

  return group_args(
    sigmaa=final_sigmaa_result.sigmaa,
    k_sol=k_sol, b_sol=b_sol,
    n_iterations=n_iterations,
    converged=converged,
    history=history)

def f_model_no_aniso_scale(fmodel):
  """ Reconstruct f_model_no_aniso_scale (the intermediate result before
  k_anisotropic is applied -- k_isotropic*(F_calc + k_mask*F_mask +
  F_part1 + F_part2), bulk solvent included -- from an mmtbx.f_model.
  manager's already-public accessors, rather than exposing a new C++/
  Python accessor on the manager itself.

  Exact (to floating-point precision: verified empirically against the
  C++ core's own f_model_no_aniso_scale_ member, mmtbx/f_model/f_model.h,
  max abs diff ~3e-14 on a real, scaled fmodel) because the C++ core
  itself computes f_model as k_anisotropic * f_model_no_aniso_scale (see
  f_model.h's core constructor), so dividing f_model() back out by
  k_anisotropic() exactly undoes that last step:
    f_model_no_aniso_scale = f_model() / k_anisotropic()

  Returns a miller.array (complex), same index set/order as fmodel.f_obs().
  """
  f_model = fmodel.f_model()
  k_aniso = fmodel.k_anisotropic()
  inv_k_aniso = 1.0 / k_aniso
  return f_model.customized_copy(data=f_model.data() * inv_k_aniso)
