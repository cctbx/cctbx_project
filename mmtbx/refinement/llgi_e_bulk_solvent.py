from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.xray import ext as xray_ext
from mmtbx import scaling
from mmtbx.refinement.llgi_sigmaa import _b_spline_design_matrix
from mmtbx.refinement.llgi_sigmaa import _spline_curvature_penalty_and_gradient
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit
import scitbx.lbfgs
from libtbx import group_args
import iotbx.phil
import mmtbx.refinement.llgi_e_dmodel_fit as llgi_e_dmodel_fit

llgi_e_bulk_solvent_params = iotbx.phil.parse("""\
  sigmaa_model = spline *d_model
    .type = choice
    .short_caption = E-scale sigmaA(resolution) functional form
    .help = "Chooses which parametrization the E-scale sigmaA(d) fit "\
            "uses:\\n"\
            "  spline   -- B-spline-over-sigmoid curve (see "\
            "e_sigmaa_target_evaluator), the original parametrization.\\n"\
            "  d_model  -- DEFAULT. Physically motivated D_model(s; theta) "\
            "parametrization from sigmaA_model_handoff.md (see doc/"\
            "llgi_target_design.md sec. 6.4 and mmtbx.refinement."\
            "llgi_e_dmodel_fit): K positive Gaussian-decay coordinate-"\
            "error terms minus one negative Gaussian bulk-solvent-"\
            "defect covariance term, fit by direct joint L-BFGS "\
            "against the full free-reflection likelihood (not a "\
            "smoothing spline). Matches the spline's R-free on 9RRL, "\
            "2G38 and 1OS8 and does not collapse at high resolution on "\
            "a poor starting model; see d_model_params for its own "\
            "sub-parameters."
  d_model_params {
    include scope mmtbx.refinement.llgi_e_dmodel_fit.llgi_e_dmodel_params
  }
  n_sigmap_nodes = 15
    .type = int
    .expert_level = 3
    .help = "Number of Chebyshev nodes (and fit terms) used for the " \
            "SigmaP(d*^2) smoothed intensity trend -- see " \
            "mmtbx.refinement.llgi_e_bulk_solvent.build_sigma_p. Matches " \
            "mmtbx.scaling.absolute_scaling.kernel_normalisation's " \
            "n_bins/n_term defaults (23/13); kept smaller and equal here " \
            "since this is recomputed for every sigmaA fit (see " \
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
            "sigmaA(resolution) curve fit against R-free data " \
            "(design note sec. 5). Same knot-" \
            "placement convention as the F-scale llgi_sigmaa module."
  spline_degree = 3
    .type = int
    .expert_level = 3
    .help = "Degree of the clamped B-spline basis used for sigmaA(d*^2)."
  sigmaa_max_iterations = 100
    .type = int
    .expert_level = 3
    .help = "Max LBFGS iterations for each sigmaA spline refit."
  sigmaa_curvature_weight = 0.02
    .type = float
    .short_caption = SigmaA spline curvature restraint weight (E-scale)
    .help = "Weight of a light restraint on the second difference of the "\
            "E-scale sigmaA B-spline's raw (pre-sigmoid) coefficients, "\
            "penalising curvature of the fitted curve -- see "\
            "mmtbx.refinement.llgi_sigmaa."\
            "_spline_curvature_penalty_and_gradient. 0 disables it."
""", process_includes=True)

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
  vanishingly unlikely), this is called for every sigmaA fit, so it adds
  one guard the original does not
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
  given -- see below). Recomputed for every sigmaA fit, since
  f_model_no_aniso_scale changes with the model and bulk solvent.

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

  Computed directly from the nacelle FEFF/RESN columns already loaded
  via phenix.refinement.llgi_data.get_llgi_data.

  Returns a flex.double, one Eeff value per input reflection.
  """
  assert resn.size() == feff.size()
  return feff / resn

def _sigmoid(z, lower=0.01, upper=0.99):
  """ Bounded sigmoid reparameterisation, matching
  mmtbx.scaling.sigmaa_estimation.sigmaa_point_estimator's convention
  (same lower/upper bounds), applied pointwise to the sigmaA B-spline
  curve. Returns (value, d(value)/d(z)), both numpy arrays.
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

def bss_k_sol_b_sol(fmodel, k_sol_default=0.35, b_sol_default=46.0):
  """ Scalar (k_sol, b_sol) summary of bss's own per-reflection k_mask()
  fit, via fmodel.k_sol_b_sol_from_k_mask() (a low-resolution Gaussian
  start refined by a local grid search, clipped to [0, 0.6]/[0, 150]).
  bss's k_mask is not in general an exact k_sol*exp(-b_sol*ss) curve, so
  this is a point estimate for logging and for D_model's B_defect anchor,
  not what Emodel is built from. Falls back to the defaults if bss has
  no estimate. Returns (k_sol, b_sol) as plain floats.
  """
  k_masks = fmodel.k_masks()
  assert len(k_masks) == 1, (
    "bss_k_sol_b_sol: only a single bulk-solvent mask shell is "
    "supported (design note sec. 9); got %d." % len(k_masks))
  k_sol, b_sol = fmodel.k_sol_b_sol_from_k_mask()
  if(k_sol is None or b_sol is None):
    return k_sol_default, b_sol_default
  return float(k_sol), float(b_sol)


class e_sigmaa_target_evaluator(object):
  """ scitbx.lbfgs target evaluator optimising the B-spline coefficients
  of the E-scale sigmaA(resolution) curve against the E-scale LLGI
  target, summed over the R-free/test set only (design note sec. 5),
  with Emodel (hence the bulk-solvent model) held fixed. The curve is a
  B-spline over a bounded sigmoid (_sigmoid), with the light curvature
  restraint mmtbx.refinement.llgi_sigmaa._spline_curvature_penalty_and_
  gradient; the target is ext.llgi_e_sigmaa_target_and_gradients.
  """

  def __init__(self,
        e_eff, test_selection, e_model, dobs, centric_flags,
        sigmaa_design, n_sigmaa_coeffs, max_iterations=100,
        curvature_weight=0.0, spline_degree=3,
      hybrid=None):
    self.hybrid = hybrid
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
      centric_flags=self.centric_flags, hybrid=self.hybrid)
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
      curvature_weight=0.0,
      hybrid=None):
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
  if(hybrid is not None):
    # exact wherever possible: a sigmaA-dependent switch would make the
    # fitted objective discontinuous (see llgi_exact.h class hybrid)
    hybrid = hybrid.with_rice_kappa(0.0)
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
    spline_degree=spline_degree, hybrid=hybrid)
  sigmaa = evaluator.sigmaa()
  final_result = xray_ext.llgi_e_sigmaa_target_and_gradients(
    e_eff=e_eff, selection=r_free_flags, e_model=e_model, dobs=dobs,
    sigmaa=sigmaa, centric_flags=centric_flags, hybrid=hybrid)
  return group_args(
    sigmaa=sigmaa, target=final_result.target(), x_range=x_range,
    evaluate_at=evaluator.evaluate_at)

def _estimate_sigmaa(e_eff, r_free_flags, e_model, dobs, centric_flags,
      d_star_sq, params, b_sol_anchor=None,
      hybrid=None):
  """ Dispatch to either the physically-motivated D_model(s; theta)
  (params.sigmaa_model=="d_model", the default -- doc/
  llgi_target_design.md sec. 6.4) or the B-spline
  (params.sigmaa_model=="spline") sigmaA(resolution) fit.

  b_sol_anchor: only used by the d_model path, where it fixes B_defect
  (see mmtbx.refinement.llgi_e_dmodel_fit.d_model_target_evaluator);
  ignored by the spline path. Callers pass bss's B_sol point estimate
  (bss_k_sol_b_sol), or None to fit B_defect instead.

  Returns a group_args with the same shape either path returns
  (.sigmaa, .target, .x_range, .evaluate_at) -- callers do not need to
  special-case which model actually ran.
  """
  if(hybrid is not None):
    # exact wherever possible: a sigmaA-dependent switch would make the
    # fitted objective discontinuous (see llgi_exact.h class hybrid)
    hybrid = hybrid.with_rice_kappa(0.0)
  if(params.sigmaa_model == "d_model"):
    dp = params.d_model_params
    return llgi_e_dmodel_fit.estimate_d_model_sigmaa(
      e_eff=e_eff, r_free_flags=r_free_flags, e_model=e_model, dobs=dobs,
      centric_flags=centric_flags, d_star_sq=d_star_sq,
      n_gaussian_terms=dp.n_gaussian_terms,
      max_iterations=dp.max_iterations,
      b_sol_anchor=b_sol_anchor,
      include_constant_term=dp.include_constant_term, hybrid=hybrid)
  return estimate_e_sigmaa(
    e_eff=e_eff, r_free_flags=r_free_flags, e_model=e_model, dobs=dobs,
    centric_flags=centric_flags, d_star_sq=d_star_sq,
    n_coeffs=params.n_sigmaa_coeffs, spline_degree=params.spline_degree,
    max_iterations=params.sigmaa_max_iterations,
    curvature_weight=params.sigmaa_curvature_weight, hybrid=hybrid)

def estimate_e_sigmaa_for_fmodel(fmodel, dobs, feff, resn, params=None):
  """ Fit the E-scale sigmaA(resolution) curve against fmodel's current
  model, with the bulk solvent exactly as bss's least-squares fit left it
  (fmodel.k_masks()[0], via f_model_no_aniso_scale; fmodel is never
  modified). Fitting bulk solvent against the E-scale LLGI instead was
  tried and dropped (doc/llgi_target_design.md sec. 0.4): it was worse
  than LS on 1OS8 and 2G38 and could pin B_sol at its bound.

  fmodel: an mmtbx.f_model.manager, already scaled (e.g. straight after
  update_all_scales()/bss).

  dobs, feff, resn: nacelle DOBS/FEFF/RESN, already matched to
  fmodel.f_obs()'s current index set.

  params: extracted llgi_e_bulk_solvent_params phil, or None for
  defaults.

  Returns a group_args with .sigmaa (flex.double, every reflection),
  .target (the fit's final test-set target), .x_range/.evaluate_at (to
  evaluate the fitted curve at other resolutions, e.g. for the fill-
  missing map step) and .k_sol/.b_sol (bss_k_sol_b_sol's point estimate
  of bss's bulk solvent, for logging; b_sol also anchors D_model's
  B_defect).
  """
  if(params is None):
    params = llgi_e_bulk_solvent_params.extract()
  f_obs = fmodel.f_obs()
  r_free_flags = fmodel.r_free_flags().data()
  centric_flags = f_obs.centric_flags().data()
  epsilons = f_obs.epsilons().data().as_double()
  d_star_sq = f_obs.d_star_sq().data()

  e_eff = build_e_eff(feff, resn)
  import mmtbx.refinement.llgi_hybrid as llgi_hybrid
  hybrid = llgi_hybrid.get_hybrid(fmodel.llgi_data())
  fmnas = f_model_no_aniso_scale(fmodel).data()
  e_model = build_e_model(
    fmnas, epsilons, d_star_sq,
    n_sigmap_nodes=params.n_sigmap_nodes,
    auto_kernel_number=params.auto_kernel_number).e_model
  k_sol, b_sol = bss_k_sol_b_sol(fmodel)
  sigmaa_result = _estimate_sigmaa(
    e_eff=e_eff, r_free_flags=r_free_flags,
    e_model=flex.abs(e_model), dobs=dobs,
    centric_flags=centric_flags, d_star_sq=d_star_sq,
    params=params, b_sol_anchor=b_sol, hybrid=hybrid)
  return group_args(
    sigmaa=sigmaa_result.sigmaa,
    target=sigmaa_result.target,
    x_range=sigmaa_result.x_range,
    evaluate_at=sigmaa_result.evaluate_at,
    k_sol=k_sol, b_sol=b_sol)

def estimate_sigmaa_e_then_scatfrac_f(
      fmodel, dobs, feff, resn, e_params=None, scatfrac_params=None):
  """ Break the F-scale sigmaA/ScatFrac non-identifiability by fitting
  sigmaA FIRST against the E-scale LLGI target, then ScatFrac against the
  F-scale LLGI target with sigmaA fixed.

  Why this ordering avoids the degeneracy without needing an empirical
  estimator: the F-scale target only sees sigmaA and ScatFrac through the
  combination D = Dobs*sigmaA/sqrt(ScatFrac), so a joint fit of both
  against it is degenerate. The E-scale target (llgi_e.h/estimate_e_
  sigmaa) has NO ScatFrac term at all -- Emodel is normalised by
  sqrt(EPS*SigmaP), not by ScatFrac -- so fitting sigmaA there first is
  not a degenerate problem, and the F-scale ScatFrac fit that follows has
  only one free curve left.

  Step 1 (sigmaA, E-scale, R-free only): delegates to
  estimate_e_sigmaa_for_fmodel, with bulk solvent as bss's LS fit left
  it.

  Step 2 (ScatFrac, F-scale, working set = all minus R-free): delegates
  to mmtbx.refinement.llgi_sigmaa.estimate_llgi_scatfrac_likelihood,
  fixing sigmaA at Step 1's result. Uses fmodel.f_model() (bulk-solvent-
  and scale-corrected), matching update_llgi_sigmaa_scatfrac's own
  convention: raw f_calc() lacks the k_isotropic correction Feff/Resn/
  f_model() already reflect. ScatFrac is NOT bounded above by 1: Feff
  need not be on absolute scale (its scaling
  assumes 50% solvent content by default), so ScatFrac can genuinely
  exceed 1 -- see llgi_scatfrac_b_factor_target_evaluator's docstring.

  fmodel: an mmtbx.f_model.manager, already scaled (e.g. straight after
  update_all_scales()/bss).

  dobs, feff, resn: nacelle DOBS/FEFF/RESN columns, already matched to
  fmodel.f_obs()'s CURRENT index set.

  e_params: extracted llgi_e_bulk_solvent_params phil, or None for
  defaults (passed straight through to estimate_e_sigmaa_fixed_bulk_
  solvent for Step 1).

  scatfrac_params: extracted llgi_sigmaa_scatfrac_params phil, or None
  for defaults (passed straight through to estimate_llgi_scatfrac_
  likelihood for Step 2).

  Returns a group_args with .sigmaa, .scatfrac (flex.double, one value
  per input reflection, evaluated at every reflection), .target (the
  Step 2 ScatFrac fit's final LLGI target value on the working set, for
  diagnostics/logging), and .scatfrac_inf/.b_scatfrac (floats,
  forwarded from estimate_llgi_scatfrac_likelihood's own result for logging/
  diagnostics, see llgi_scatfrac_b_factor_target_evaluator's docstring
  for why these two values should not be over-interpreted on their own).
  """
  import mmtbx.refinement.llgi_sigmaa as llgi_sigmaa
  if(e_params is None):
    e_params = llgi_e_bulk_solvent_params.extract()
  if(scatfrac_params is None):
    scatfrac_params = llgi_sigmaa.llgi_sigmaa_scatfrac_params.extract()

  sigmaa_result = estimate_e_sigmaa_for_fmodel(
    fmodel, dobs=dobs, feff=feff, resn=resn, params=e_params)
  sigmaa = sigmaa_result.sigmaa

  f_obs = fmodel.f_obs()
  r_free_flags = fmodel.r_free_flags().data()
  working_selection = ~r_free_flags
  llgi_data = fmodel.llgi_data()
  teps = llgi_data.teps.data()
  import mmtbx.refinement.llgi_hybrid as llgi_hybrid
  hybrid = llgi_hybrid.get_hybrid(llgi_data)
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
    params=scatfrac_params, hybrid=hybrid)

  return group_args(
    sigmaa=sigmaa, scatfrac=scatfrac_result.scatfrac,
    target=scatfrac_result.target,
    scatfrac_inf=getattr(scatfrac_result, "scatfrac_inf", None),
    b_scatfrac=getattr(scatfrac_result, "b_scatfrac", None))

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
