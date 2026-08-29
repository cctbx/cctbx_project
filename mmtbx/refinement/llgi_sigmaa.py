from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.xray import ext
import scitbx.lbfgs
from libtbx import group_args
import iotbx.phil

llgi_sigmaa_scatfrac_params = iotbx.phil.parse("""\
  n_sigmaa_coeffs = 8
    .type = int
    .short_caption = Number of sigmaA spline coefficients
    .help = "Number of B-spline coefficients for the sigmaA(resolution) " \
            "curve. Knots are placed evenly in d*^2."
  n_scatfrac_coeffs = 8
    .type = int
    .short_caption = Number of ScatFrac spline coefficients
    .help = "Number of B-spline coefficients for the ScatFrac(resolution) " \
            "curve. Same d*^2 knot placement convention as sigmaA, " \
            "independent coefficients."
  spline_degree = 3
    .type = int
    .expert_level = 3
    .help = "Degree of the clamped B-spline basis used for both curves."
  n_scatfrac_bins = 20
    .type = int
    .expert_level = 3
    .help = "Number of d*^2 bins ScatFrac's per-reflection Fcalc^2/" \
            "(Teps*Resn^2) ratio is averaged over before spline-fitting " \
            "in log space, weighted by each bin's effective reflection " \
            "count (n_acentric + n_centric/2). See " \
            "mmtbx.refinement.llgi_sigmaa.estimate_llgi_scatfrac."
  max_iterations = 100
    .type = int
    .expert_level = 3
""")

def _b_spline_design_matrix(x, n_coeffs, degree, x_range=None):
  """ Clamped B-spline design matrix B[i,k] = B_k(x[i]), for a basis with
  interior knots evenly spaced over x_range (mapped internally to [0,1]
  and back, so this works for any x range including d*^2). Returns a
  numpy array of shape (len(x), n_coeffs). Since a curve z(x) = design @
  coeffs is linear in coeffs, this matrix needs to be built only once per
  macrocycle (x -- the resolution metric per reflection -- does not
  change during a sigmaA/ScatFrac fit) and is reused every LBFGS
  iteration (for sigmaA) or in a single least-squares solve (for
  ScatFrac); see doc/llgi_target_design.md sec. 5.2.

  x_range: (x_min, x_max) to normalise against. If None, derived from
  x itself (the historical default, fine when this function is called
  once against the full reflection set). MUST be passed explicitly and
  consistently when this function is called more than once for the same
  fit against different x arrays (e.g. once against per-bin centres to
  solve for spline coefficients, once against the full per-reflection
  array to evaluate the fitted curve) -- otherwise the two calls
  normalise x differently (bin centres never span the full data's
  min/max) and the resulting coefficients get silently misapplied when
  evaluated on the second array. Caught while adding per-bin fitting to
  estimate_llgi_scatfrac, which calls this function twice for exactly
  that reason.
  """
  import numpy as np
  from scipy.interpolate import BSpline
  x = np.asarray(x, dtype=float)
  if(x_range is None):
    x_min = float(x.min())
    x_max = float(x.max())
  else:
    x_min, x_max = x_range
  if(x_max <= x_min):
    # Degenerate resolution range (e.g. a single reflection, or all
    # reflections at identical resolution): fall back to a constant
    # basis (every reflection maps to the same single coefficient) so
    # this does not crash; the caller ends up fitting one overall value.
    x_max = x_min + 1.0
  x_norm = (x - x_min) / (x_max - x_min)
  n_knots = n_coeffs + degree + 1
  n_interior = n_knots - 2 * degree
  if(n_interior < 2):
    raise RuntimeError(
      "n_coeffs=%d is too small for spline_degree=%d (need at least "
      "%d coefficients)." % (n_coeffs, degree, 2 * degree - degree + 1))
  interior = np.linspace(0.0, 1.0, n_interior)
  knots = np.concatenate([[0.0] * degree, interior, [1.0] * degree])
  design = BSpline.design_matrix(
    x_norm, knots, degree, extrapolate=False).toarray()
  return design

def _sigmoid(z, lower=0.01, upper=0.99):
  """ Bounded sigmoid reparameterisation, matching
  mmtbx.scaling.sigmaa_estimation.sigmaa_point_estimator's convention
  exactly (same lower/upper bounds), applied pointwise to a B-spline
  curve evaluated at z rather than to a single scalar LBFGS parameter.
  Returns (value, d(value)/d(z)), both as numpy arrays.
  """
  import numpy as np
  z = np.asarray(z, dtype=float)
  exp_neg_z = np.exp(-z)
  value = lower + (upper - lower) / (1.0 + exp_neg_z)
  dvalue_dz = (upper - lower) * exp_neg_z / (1.0 + exp_neg_z) ** 2
  return value, dvalue_dz

def estimate_llgi_scatfrac(
      f_calc,
      teps,
      resn,
      d_star_sq,
      centric_flags,
      scale_factor=1.0,
      n_coeffs=8,
      spline_degree=3,
      n_bins=20):
  """ Fit ScatFrac(resolution) as a B-spline curve by weighted least-
  squares to per-bin RATIO-OF-SUMS estimates (not a mean of per-
  reflection ratios -- see "Binning, point estimate, and weighting"
  below), using sum(scale_factor*Fcalc_i)^2)/sum(Teps_i*Resn_i^2) within
  each bin, whose expectation, by construction of Teps/Resn (see llgi.h
  and doc/llgi_target_design.md sec. 4.2), is ScatFrac -- the fraction of
  total scattering accounted for by the model as a function of
  resolution. This is a direct empirical calculation, NOT an LLGI-
  likelihood fit: sigmaA and ScatFrac are not jointly identifiable from
  the LLGI target alone (D = Dobs*sigmaA/sqrt(ScatFrac) is the only
  combined quantity the target sees, so a joint LBFGS fit of both curves
  has a degenerate ridge along which target value is unchanged --
  confirmed empirically before this design was adopted; see doc/
  llgi_target_design.md sec. 5.2). Computing ScatFrac this way instead
  breaks that degeneracy: sigmaA (see estimate_llgi_sigmaa) is then fit
  against LLGI with ScatFrac already fixed.

  Binning, point estimate, and weighting (the point estimate itself was
  corrected after both an earlier unweighted per-reflection fit, AND a
  later per-bin MEAN-of-ratios fit, were each found, on real data, to
  badly overestimate ScatFrac at high resolution -- see doc/
  llgi_target_design.md sec. 5.2.2/5.2.3 for the full account of both):
  the per-reflection ratio r_i = (scale_factor*Fcalc_i)^2/(Teps_i*
  Resn_i^2) is, under a Wilson-distribution assumption for |Fcalc|^2
  within a narrow resolution range, itself Wilson/exponential(-like)
  distributed with Var(r_i) = ScatFrac^2 for acentric reflections and
  2*ScatFrac^2 for centric (a standard factor-of-2 relationship) -- i.e.
  large individual values are not outliers, they are the expected heavy-
  tailed shape of the distribution. Averaging r_i directly (mean(X/Y))
  lets a single extreme |Fcalc_i| dominate the bin average on its own,
  since squaring happens before any averaging; this was observed
  directly on real data -- a single reflection whose |Fcalc| swung ~10x
  between two macrocycles of ordinary refinement inflated an entire
  ~1600-reflection bin's mean ratio by ~10x on its own (bin median was
  unaffected). The fix is to sum numerator and denominator SEPARATELY
  across the bin first, then divide (ratio(sum(X), sum(Y)), equivalently
  the ratio of the bin means) -- mathematically distinct from mean(X/Y)
  for a heavy-tailed X, and far more robust to exactly this kind of
  single-reflection swing, since the sum of many Fcalc^2 values averages
  out before the division happens. This function: (1) bins reflections
  into n_bins equal-population-ish bins by d_star_sq (numpy.array_split
  on the sorted metric), (2) computes each bin's ratio-of-sums estimate
  r_bar and an effective sample size n_eff = n_acentric + n_centric/2
  (correcting for centric reflections' 2x variance), (3) fits the
  B-spline in LOG space (z(x) = ln ScatFrac(x)) rather than to ScatFrac
  directly -- by the delta method, Var(ln r_bar) ~= 1/n_eff (re-derived
  and verified numerically for the ratio-of-sums estimator specifically,
  confirming the 1/n scaling carries over unchanged from the original
  single-reflection-ratio derivation -- only the point estimate itself
  needed to change, not the weighting), so weighting by n_eff directly
  in log-space is both simpler and better-motivated than carrying
  Var(r_bar) = ScatFrac^2/n_eff (circular, since ScatFrac is the
  unknown) through an untransformed fit. Fitting in log-space also
  guarantees ScatFrac > 0 automatically, a loose end an earlier version
  needed an ad hoc clip for.

  scale_factor: the same overall scale factor k used elsewhere in the
  LLGI target (llgi.h's `k`, typically manager.scale_ml_wrapper()) --
  Fcalc from an xray_structure/fmodel manager is not guaranteed to
  already be on the same absolute scale as Feff/Resn (which derive from
  the experimental data), so it must be rescaled by k for this ratio to
  mean what "fraction of total scattering, ~1 for a complete model"
  requires. Omitting this (an earlier draft of this function did) can
  give wildly implausible ScatFrac values whenever Fcalc's absolute scale
  differs from Feff/Resn's, as was caught testing against a real (if
  crudely scaled) mmtbx.f_model.manager rather than only hand-constructed
  unit-scale test arrays.

  Uses the FULL reflection set (not R-free/test-set-only), since this is
  an empirical intensity ratio, not a likelihood optimisation at risk of
  overfitting to the set it is validated against.

  Returns a flex.double, one ScatFrac value per input reflection
  (evaluated at that reflection's own d_star_sq).
  """
  import numpy as np
  teps_np = teps.as_numpy_array()
  resn_np = resn.as_numpy_array()
  fc_abs_sq = ((flex.abs(f_calc) * scale_factor) ** 2).as_numpy_array()
  denom = teps_np * resn_np ** 2
  d_star_sq_np = d_star_sq.as_numpy_array()
  is_centric = centric_flags.as_numpy_array()

  order = np.argsort(d_star_sq_np)
  n_bins_eff = max(1, min(n_bins, len(order)))
  bin_indices = np.array_split(order, n_bins_eff)
  bin_x = []
  bin_log_ratio = []
  bin_weight = []
  for idx in bin_indices:
    if(len(idx) == 0):
      continue
    # Ratio of bin SUMS (sum(Fcalc^2) / sum(Teps*Resn^2)), NOT the mean of
    # the n per-reflection ratios Fcalc_i^2/(Teps_i*Resn_i^2). These are
    # mathematically different quantities for a heavy-tailed numerator
    # (mean(X/Y) != mean(X)/mean(Y) in general), and the per-reflection-
    # ratio-mean version was found, on real data, to be badly distorted
    # by a small number of individual reflections with large |Fcalc| --
    # squaring in the denominator (small Resn at high resolution)
    # inflates that single reflection's own ratio term directly, before
    # any averaging happens, so it dominates the bin mean even at n~1600
    # (traced to a specific reflection whose |Fcalc| happened to swing
    # ~10x between two macrocycles of real refinement -- a real, if
    # extreme, per-reflection Fcalc change, not a bug in Fcalc/Resn/Teps
    # themselves; the bug was averaging the ratio rather than the ratio
    # of averages). The sum-ratio construction lets Fcalc^2 values
    # average out *before* dividing by the (roughly constant within a
    # narrow bin) denominator, which is far more robust to exactly this
    # kind of outlier. Re-derived and numerically verified (Monte Carlo)
    # that Var(sum-ratio) = ScatFrac^2/n, the same 1/n scaling as the
    # original (single-reflection-ratio) derivation, so the n_eff
    # weighting and log-space delta-method argument below are unaffected
    # by this fix -- only the point estimate itself needed to change.
    r_bar = fc_abs_sq[idx].sum() / denom[idx].sum()
    if(r_bar <= 0):
      continue
    n_acentric = np.count_nonzero(~is_centric[idx])
    n_centric = np.count_nonzero(is_centric[idx])
    n_eff = n_acentric + 0.5 * n_centric
    if(n_eff <= 0):
      continue
    bin_x.append(d_star_sq_np[idx].mean())
    bin_log_ratio.append(np.log(r_bar))
    bin_weight.append(n_eff)
  bin_x = np.array(bin_x)
  bin_log_ratio = np.array(bin_log_ratio)
  bin_weight = np.array(bin_weight)

  # Both design matrices below MUST share the same normalisation range
  # (the full dataset's d*^2 range, not the narrower range spanned by
  # bin centres) -- see _b_spline_design_matrix's x_range docstring.
  x_range = (float(d_star_sq_np.min()), float(d_star_sq_np.max()))
  design_bins = _b_spline_design_matrix(
    bin_x, n_coeffs, spline_degree, x_range=x_range)
  sqrt_w = np.sqrt(bin_weight)
  coeffs, _residuals, _rank, _sv = np.linalg.lstsq(
    design_bins * sqrt_w[:, None], bin_log_ratio * sqrt_w, rcond=None)

  design_all = _b_spline_design_matrix(
    d_star_sq_np, n_coeffs, spline_degree, x_range=x_range)
  fitted_log = design_all.dot(coeffs)
  fitted = np.exp(fitted_log)
  return flex.double(fitted)

class llgi_sigmaa_target_evaluator(object):
  """ scitbx.lbfgs target evaluator optimising the B-spline coefficients
  of sigmaA(resolution) against the LLGI target, summed over the R-free/
  test set only (see doc/llgi_target_design.md sec. 5.2), with ScatFrac
  held fixed (see estimate_llgi_scatfrac -- sigmaA and ScatFrac are not
  jointly identifiable from LLGI alone, so ScatFrac must already be
  determined before this runs). x is the (unconstrained) B-spline
  coefficient vector; sigmaA itself is bounded to (0.01, 0.99) via the
  sigmoid reparameterisation in _sigmoid, applied after spline
  evaluation, matching sigmaa_estimation.py's convention.

  Follows scitbx.lbfgs conventions directly (see
  mmtbx/scaling/sigmaa_estimation.py's sigmaa_point_estimator): x is a
  flex.double parameter vector mutated in place by scitbx.lbfgs.run,
  compute_functional_and_gradients() returns (f, g). The underlying
  cctbx.xray.targets.llgi target (llgi.h) is already in this codebase's
  minimize-me convention (see llgi.h's target_one_h docstring), so no
  extra sign flip is applied here beyond the one already baked into the
  C++ target.
  """

  def __init__(self,
        f_eff,
        selection,
        f_calc,
        dobs,
        scatfrac,
        teps,
        resn,
        centric_flags,
        scale_factor,
        sigmaa_design,
        n_sigmaa_coeffs,
        max_iterations=100):
    self.f_eff = f_eff
    self.selection = selection
    self.f_calc = f_calc
    self.dobs = dobs
    self.scatfrac = scatfrac
    self.teps = teps
    self.resn = resn
    self.centric_flags = centric_flags
    self.scale_factor = scale_factor
    self.sigmaa_design = sigmaa_design  # numpy array, (n_refl, n_sigmaa_coeffs)
    self.n_sigmaa_coeffs = n_sigmaa_coeffs
    # Unconstrained starting point: z=0 maps (via the sigmoid) to
    # sigmaA=0.5, a neutral starting guess.
    self.x = flex.double(n_sigmaa_coeffs, 0.0)
    term_parameters = scitbx.lbfgs.termination_parameters(
      max_iterations=max_iterations)
    # As in sigmaa_point_estimator: the sigmoid reparameterisation makes
    # the line search behave poorly right at the (0.01, 0.99) edges.
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
    result = ext.llgi_sigmaa_scatfrac_target_and_gradients(
      f_eff=self.f_eff,
      selection=self.selection,
      f_calc=self.f_calc,
      dobs=self.dobs,
      sigmaa=flex.double(sigmaa),
      scatfrac=self.scatfrac,
      scale_factor=self.scale_factor,
      teps=self.teps,
      resn=self.resn,
      centric_flags=self.centric_flags)
    f = result.target()
    d_target_by_dsigmaa = np.array(result.d_target_by_dsigmaa())
    # Chain rule: d(target)/d(coeffs) = design.T @ (d_target * dsigmoid_dz),
    # since z_i = design[i,:] . coeffs is linear in coeffs.
    g = self.sigmaa_design.T.dot(d_target_by_dsigmaa * dsigmaa_dz)
    return f, flex.double(g)

  def sigmaa(self):
    """ Return the fitted sigmaA array at the current (final, once the
    minimizer has run) coefficient values. """
    sigmaa, _ = self._current_sigmaa()
    return flex.double(sigmaa)

def estimate_llgi_sigmaa(
      f_eff,
      r_free_flags,
      f_calc,
      dobs,
      scatfrac,
      teps,
      resn,
      centric_flags,
      d_star_sq,
      scale_factor=1.0,
      params=None):
  """ Fit sigmaA(resolution) as a B-spline curve against the LLGI target,
  restricted to the R-free/test set, with ScatFrac already fixed (see
  estimate_llgi_scatfrac -- must be called first; sigmaA and ScatFrac are
  not jointly identifiable from LLGI alone). Evaluates the fitted curve at
  every reflection (not just the test set) so the result can be used
  directly by the main llgi refinement target.

  params: extracted llgi_sigmaa_scatfrac_params phil, or None for
  defaults (only n_sigmaa_coeffs, spline_degree, max_iterations are used
  here; n_scatfrac_coeffs is used by estimate_llgi_scatfrac).

  Returns a group_args with .sigmaa (flex.double, one value per input
  reflection) and .target (final fitted LLGI target value on the test
  set, for diagnostics/logging).
  """
  if(params is None):
    params = llgi_sigmaa_scatfrac_params.extract()
  n_refl = f_eff.size()
  assert r_free_flags.size() == n_refl
  assert f_calc.size() == n_refl
  assert dobs.size() == n_refl
  assert scatfrac.size() == n_refl
  assert teps.size() == n_refl
  assert resn.size() == n_refl
  assert centric_flags.size() == n_refl
  assert d_star_sq.size() == n_refl
  n_test = r_free_flags.count(True)
  if(n_test == 0):
    raise RuntimeError(
      "No R-free/test-set reflections available for the LLGI sigmaA fit.")
  sigmaa_design = _b_spline_design_matrix(
    d_star_sq.as_numpy_array(), params.n_sigmaa_coeffs, params.spline_degree)
  evaluator = llgi_sigmaa_target_evaluator(
    f_eff=f_eff,
    selection=r_free_flags,
    f_calc=f_calc,
    dobs=dobs,
    scatfrac=scatfrac,
    teps=teps,
    resn=resn,
    centric_flags=centric_flags,
    scale_factor=scale_factor,
    sigmaa_design=sigmaa_design,
    n_sigmaa_coeffs=params.n_sigmaa_coeffs,
    max_iterations=params.max_iterations)
  sigmaa = evaluator.sigmaa()
  final_result = ext.llgi_sigmaa_scatfrac_target_and_gradients(
    f_eff=f_eff, selection=r_free_flags, f_calc=f_calc, dobs=dobs,
    sigmaa=sigmaa, scatfrac=scatfrac, scale_factor=scale_factor,
    teps=teps, resn=resn, centric_flags=centric_flags)
  return group_args(sigmaa=sigmaa, target=final_result.target())

def estimate_llgi_sigmaa_scatfrac(
      f_eff,
      r_free_flags,
      f_calc,
      dobs,
      teps,
      resn,
      centric_flags,
      d_star_sq,
      scale_factor=1.0,
      params=None):
  """ Convenience wrapper: computes ScatFrac(resolution) empirically
  (estimate_llgi_scatfrac, full reflection set) then fits sigmaA(
  resolution) against LLGI with ScatFrac held fixed (estimate_llgi_sigmaa,
  R-free/test set only). See those two functions' docstrings for why this
  two-step (not jointly-optimised) design was adopted.

  Returns a group_args with .sigmaa, .scatfrac (flex.double, one value
  per input reflection) and .target (final fitted LLGI target value on
  the test set, for diagnostics/logging).
  """
  if(params is None):
    params = llgi_sigmaa_scatfrac_params.extract()
  scatfrac = estimate_llgi_scatfrac(
    f_calc=f_calc, teps=teps, resn=resn, d_star_sq=d_star_sq,
    centric_flags=centric_flags,
    scale_factor=scale_factor,
    n_coeffs=params.n_scatfrac_coeffs, spline_degree=params.spline_degree,
    n_bins=params.n_scatfrac_bins)
  sigmaa_result = estimate_llgi_sigmaa(
    f_eff=f_eff, r_free_flags=r_free_flags, f_calc=f_calc, dobs=dobs,
    scatfrac=scatfrac, teps=teps, resn=resn, centric_flags=centric_flags,
    d_star_sq=d_star_sq, scale_factor=scale_factor, params=params)
  return group_args(
    sigmaa=sigmaa_result.sigmaa, scatfrac=scatfrac,
    target=sigmaa_result.target)
