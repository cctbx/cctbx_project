from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.xray import ext
import scitbx.lbfgs
from libtbx import group_args
import iotbx.phil

llgi_sigmaa_scatfrac_params = iotbx.phil.parse("""\
  max_iterations = 100
    .type = int
    .expert_level = 3
  scatfrac_b_factor_restraint_sigma = 10.0
    .type = float
    .short_caption = ScatFrac B-factor restraint sigma (Angstrom^2)
    .help = "Weak quadratic restraint on B_scatfrac toward 0, in the "\
            "ScatFrac_inf*exp(-B_scatfrac*ss) ScatFrac curve. "\
            "Motivation: sigmaA and "\
            "ScatFrac are not independently identifiable -- the F-scale "\
            "LLGI target constrains only their combination -- so B_scatfrac "\
            "and sigmaA can trade off against each other over the observed "\
            "resolution range with no likelihood cost (see llgi_scatfrac_"\
            "b_factor_target_evaluator's docstring on the ScatFrac_inf/"\
            "B_scatfrac extrapolation ambiguity, the same phenomenon one "\
            "level up). This degeneracy is not perfectly free, though: "\
            "sigmaA is bounded to [0, 1], so a B_scatfrac far from 0 that "\
            "can only be compensated by pushing sigmaA toward or past its "\
            "upper bound is not actually a cost-free direction, just one "\
            "the unrestrained fit may wander into anyway before that "\
            "bound is felt (e.g. mid-refinement, far from convergence). "\
            "Restraining B_scatfrac directly -- rather than sigmaA, which "\
            "is already bounded and fit by a separate, earlier E-scale "\
            "step -- damps that wandering at its source without touching "\
            "the sigmaA fit at all. Penalty is 0.5*(B_scatfrac/sigma)^2 "\
            "(a standard restrained-B-factor form, sigma in the same "\
            "Angstrom^2 units as B_scatfrac itself): weak enough to be "\
            "negligible next to the LLGI likelihood whenever the data "\
            "themselves prefer a real trend, but enough to pull B_scatfrac "\
            "back toward the flat (scalar-like) fit when the likelihood is "\
            "close to degenerate over the observed range. See "\
            "_b_factor_restraint_penalty_and_gradient for the mechanism. "\
            "0 disables it."
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

def _spline_curvature_penalty_and_gradient(coeffs, weight):
  """ Light restraint discouraging curvature in a B-spline curve's raw
  (pre-sigmoid, "z-space") coefficients -- added to a sigmaA(resolution)
  LLGI fit's target/gradient to stop the curve collapsing sharply toward
  its lower bound in resolution shells where the R-free test set is too
  sparse for the LLGI likelihood alone to constrain the fit (observed on
  real data: sigmaA(d) plunging to the sigmoid floor over the last couple
  of resolution shells, well before the reflections actually run out,
  then recovering somewhat on the next macrocycle -- not physically
  expected, since fit quality should vary smoothly with resolution).

  Penalises the discrete second difference of the coefficient vector,
  R(c) = weight * sum_i (c[i-1] - 2*c[i] + c[i+1])^2, a standard roughness
  penalty approximating the curvature (second derivative) of z(x) =
  design(x) . c in the spline's normalised-[0,1] domain (see
  _b_spline_design_matrix). Interior knots are placed evenly in that
  domain, so a constant knot spacing is implicit in treating the second
  difference of c as a curvature proxy for z itself.

  Deliberately restrains z (the pre-sigmoid quantity), not log(sigmaA)
  directly: near the sigmoid's lower bound `lower` (where the collapse
  this is meant to fix actually happens), sigmaA - lower ~= (upper-lower)
  *exp(z) for very negative z, so log(sigmaA - lower) ~= const + z there
  -- i.e. z itself is already approximately log-linear in exactly the
  regime this restraint targets, without needing to differentiate
  through the sigmoid (which would make the penalty non-quadratic in c
  and require the sigmoid's second derivative too). In the well-
  determined middle of the resolution range z is simply the natural
  unconstrained LBFGS parameter, so penalising its curvature there is
  the ordinary smoothing-spline move.

  Quadratic in c => the penalty vanishes exactly for any c that is a
  straight line (or constant) in index space (zero second difference),
  so a genuinely log-linear-like fit is entirely unaffected; it grows
  only where the fit curves, and the growth is independent of how well-
  determined that curvature is by the data -- so a small fixed weight
  has negligible relative effect where the LLGI target's own curvature
  (from many reflections) dominates, and a comparatively larger relative
  effect exactly where that likelihood curvature is weak (few/noisy
  high-resolution reflections). This gets most of the benefit of an
  adaptive (curvature- or reflection-count-weighted) penalty without
  needing one; see doc/llgi_target_design.md sec. 6 for the fuller
  discussion of alternatives considered (a hard monotonicity restraint
  was rejected: a real bulk-solvent-incomplete low-resolution rise-then-
  fall in sigmaA is physically expected and would be fought by a
  monotonicity restraint).

  coeffs: numpy array, the current (unconstrained) B-spline coefficient
  vector (i.e. target_evaluator.x as a numpy array -- NOT sigmaA itself).
  weight: penalty weight (llgi_sigmaa_scatfrac_params.sigmaa_curvature_
  weight or the E-scale equivalent); 0 (or coeffs.size() < 3) returns a
  no-op (0.0, zeros).

  Returns (penalty, d(penalty)/d(coeffs)), penalty a plain float and the
  gradient a numpy array the same shape as coeffs, both ready to be added
  directly onto compute_functional_and_gradients()'s (f, g).
  """
  import numpy as np
  n = coeffs.shape[0]
  if(weight == 0 or n < 3):
    return 0.0, np.zeros_like(coeffs)
  d2 = coeffs[:-2] - 2.0 * coeffs[1:-1] + coeffs[2:]
  penalty = weight * float(np.sum(d2 * d2))
  grad = np.zeros_like(coeffs)
  grad[:-2]  += 2.0 * weight * d2
  grad[1:-1] += -4.0 * weight * d2
  grad[2:]   += 2.0 * weight * d2
  return penalty, grad

def _b_factor_restraint_penalty_and_gradient(b_scatfrac, sigma):
  """ Weak quadratic restraint pulling B_scatfrac (llgi_scatfrac_b_factor_
  target_evaluator's slope parameter) toward 0, standard restrained-
  B-factor form: penalty = 0.5*(b_scatfrac/sigma)^2, so d(penalty)/
  d(b_scatfrac) = b_scatfrac/sigma^2. See llgi_sigmaa_scatfrac_params.
  scatfrac_b_factor_restraint_sigma's help for the motivation (damping
  the sigmaA/B_scatfrac trade-off's tendency to wander before sigmaA's
  own [0, 1] bound is reached).

  sigma: restraint width in the same Angstrom^2 units as b_scatfrac;
  <= 0 (matching 0 meaning "off") returns a no-op (0.0, 0.0).

  Returns (penalty, d(penalty)/d(b_scatfrac)), both plain floats, ready
  to be added directly onto the (f, g) pair compute_functional_and_
  gradients() returns (g's B_scatfrac component only -- ScatFrac_inf is
  not restrained by this term).
  """
  if(sigma <= 0):
    return 0.0, 0.0
  penalty = 0.5 * (b_scatfrac / sigma) ** 2
  grad = b_scatfrac / (sigma * sigma)
  return penalty, grad

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
  likelihood fit; estimate_llgi_scatfrac_likelihood uses its single-bin
  (n_bins=1) value as the starting point for the likelihood fit.

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

class llgi_scatfrac_b_factor_target_evaluator(object):
  """ scitbx.lbfgs target evaluator optimising a TWO-PARAMETER ScatFrac
  curve, ScatFrac_inf * exp(-B_scatfrac * ss) (ss = d*^2/4, the same
  sin(theta)/lambda squared convention as the bulk-solvent k_mask formula,
  so
  B_scatfrac is directly comparable to an ordinary crystallographic
  B-factor, positive B meaning ScatFrac falls off toward high resolution
  as usual), against the F-scale LLGI target, sigmaA held fixed.
  Earlier B-spline and single-value forms were dropped: a weight sweep
  of the spline's curvature restraint found no safe operating point on
  real data, and a single value cannot represent resolution dependence.

  Motivation: resolution dependence in the fraction of scattering
  accounted for is physically plausible in more than one
  direction: a partial model built from its best-ordered components
  first (as coordinate refinement of an initially poor model often
  proceeds) would have ScatFrac RISE toward high resolution (the model
  accounts for a larger share of the weak, well-ordered high-resolution
  data than of the disordered low-resolution features it has not yet
  captured) -- the opposite sign from the more familiar case of ScatFrac
  falling with resolution. A single positive-only B-factor-style
  parameter B_scatfrac (sign unconstrained here, unlike an ordinary ADP)
  covers both directions with only one more parameter than a constant.

  Design matrix is FIXED and explicit ([1, -ss], not a B-spline basis):
  z(ss) = z_inf + (-B_scatfrac/4) * d_star_sq = z_inf - B_scatfrac * ss,
  so the two coefficients are directly interpretable (z_inf = log
  ScatFrac_inf, the d*^2 -> 0 / infinite-resolution limit; the slope
  coefficient is -B_scatfrac/4) -- unlike reusing _b_spline_design_
  matrix with n_coeffs=2, degree=1, whose two coefficients are the
  spline's values at the two ENDS of the observed d*^2 range (an affine
  reparameterisation of the same line, but not directly z_inf/B_scatfrac,
  and silently redefined if the observed resolution range changes
  between macrocycles, since that spline normalises against the observed
  min/max rather than against d*^2=0).

  Log-space reparameterisation (ScatFrac = exp(z), NO upper bound of 1:
  Feff need not be on absolute scale, so ScatFrac can exceed 1). A
  straight line in log-ScatFrac-vs-ss space has no curvature to
  restrain.

  IMPORTANT: ScatFrac_inf (the d*^2 -> 0 EXTRAPOLATED intercept) and
  B_scatfrac trade off against each other over any finite resolution
  range, much like the classic Wilson-plot scale/B-factor ambiguity --
  many (ScatFrac_inf, B_scatfrac) pairs give nearly the SAME fitted curve
  over the OBSERVED range while extrapolating to very different d*^2=0
  values. Confirmed directly: a synthetic recovery test with a true
  (ScatFrac_inf, B_scatfrac) = (0.95, 25.0), observed over d*^2 in
  [0.001, 0.25], recovered (0.50, 10.5) -- individually far from the
  truth -- yet the fitted CURVE matched the true curve over that same
  observed range to 22% mean relative difference with log-log
  correlation 0.9999997. This is expected, not a bug: the two
  parameters' individual values should not be over-interpreted (e.g.
  logged/compared across macrocycles) without also checking the fitted
  curve's shape over the actual resolution range in use; only
  B_scatfrac's SIGN (falling vs. rising trend) is a robust, coarse
  property largely unaffected by this ambiguity.

  A second, related trade-off is against sigmaA rather than against
  ScatFrac_inf: sigmaA and ScatFrac are only jointly identifiable
  through the F-scale LLGI target (see mmtbx.refinement.
  llgi_e_bulk_solvent.estimate_sigmaa_e_then_scatfrac_f), so a B_scatfrac far
  from 0 can sometimes be compensated by pushing sigmaA toward its own
  upper bound of 1 with little likelihood cost -- except where that
  compensation would require sigmaA > 1, which is not available, so
  the fit is not actually free to wander that far; it is only free to
  wander until sigmaA's bound is felt, which may be well past a
  physically reasonable B_scatfrac. restraint_sigma (Angstrom^2, see
  llgi_sigmaa_scatfrac_params.scatfrac_b_factor_restraint_sigma) adds
  a weak quadratic penalty pulling B_scatfrac back toward 0 -- see
  _b_factor_restraint_penalty_and_gradient -- damping that wandering
  directly rather than relying on sigmaA's bound to arrest it.
  """

  def __init__(self,
        f_eff,
        selection,
        f_calc,
        dobs,
        sigmaa,
        teps,
        resn,
        ss,
        centric_flags,
        scale_factor,
        scatfrac_inf_start,
        b_scatfrac_start,
        restraint_sigma=0.0,
        max_iterations=100,
      hybrid=None):
    self.hybrid = hybrid
    self.f_eff = f_eff
    self.selection = selection
    self.f_calc = f_calc
    self.dobs = dobs
    self.sigmaa = sigmaa
    self.teps = teps
    self.resn = resn
    self.centric_flags = centric_flags
    self.scale_factor = scale_factor
    self.restraint_sigma = restraint_sigma
    self.n_refl = f_eff.size()
    import numpy as np
    import math
    self.ss = np.asarray(ss, dtype=float)
    assert self.ss.shape[0] == self.n_refl
    z_inf0 = math.log(float(scatfrac_inf_start))
    self.x = flex.double([z_inf0, float(b_scatfrac_start)])
    term_parameters = scitbx.lbfgs.termination_parameters(
      max_iterations=max_iterations)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=term_parameters)

  def _current_scatfrac(self):
    import numpy as np
    z_inf, b_scatfrac = self.x[0], self.x[1]
    z = z_inf - b_scatfrac * self.ss
    scatfrac = np.exp(z)
    return scatfrac  # dscatfrac/dz == scatfrac itself (per-reflection)

  def compute_functional_and_gradients(self):
    import numpy as np
    scatfrac = self._current_scatfrac()
    result = ext.llgi_sigmaa_scatfrac_target_and_gradients(
      f_eff=self.f_eff,
      selection=self.selection,
      f_calc=self.f_calc,
      dobs=self.dobs,
      sigmaa=self.sigmaa,
      scatfrac=flex.double(scatfrac),
      scale_factor=self.scale_factor,
      teps=self.teps,
      resn=self.resn,
      centric_flags=self.centric_flags, hybrid=self.hybrid)
    f = result.target()
    d_target_by_dscatfrac = np.array(result.d_target_by_dscatfrac())
    # Chain rule: d(target)/dz_inf = sum_i d_target_by_dscatfrac[i] *
    # dscatfrac_i/dz_inf, with dscatfrac_i/dz_inf = scatfrac_i (since
    # z_i = z_inf - B*ss_i, dz_i/dz_inf = 1); d(target)/dB = sum_i
    # d_target_by_dscatfrac[i] * scatfrac_i * dz_i/dB, with dz_i/dB =
    # -ss_i.
    dtarget_dscatfrac_times_scatfrac = d_target_by_dscatfrac * scatfrac
    g_z_inf = float(np.sum(dtarget_dscatfrac_times_scatfrac))
    g_b = float(np.sum(dtarget_dscatfrac_times_scatfrac * (-self.ss)))
    b_scatfrac = self.x[1]
    restraint_penalty, restraint_grad = \
      _b_factor_restraint_penalty_and_gradient(
        b_scatfrac, self.restraint_sigma)
    f += restraint_penalty
    g_b += restraint_grad
    return f, flex.double([g_z_inf, g_b])

  def scatfrac(self):
    """ Return the fitted ScatFrac curve, one value per input
    reflection. """
    return flex.double(self._current_scatfrac())

  def scatfrac_inf_and_b(self):
    """ Return (ScatFrac_inf, B_scatfrac) at the current (final, once the
    minimizer has run) parameter values -- for logging/diagnostics. """
    import math
    z_inf, b_scatfrac = self.x[0], self.x[1]
    return math.exp(z_inf), b_scatfrac

def estimate_llgi_scatfrac_likelihood(
      f_eff,
      working_selection,
      f_calc,
      dobs,
      sigmaa,
      teps,
      resn,
      centric_flags,
      d_star_sq,
      scale_factor=1.0,
      params=None,
      hybrid=None):
  """ Fit ScatFrac(resolution) = ScatFrac_inf*exp(-B_scatfrac*ss) against
  the F-scale LLGI target on the working set (all reflections minus
  R-free), with sigmaA already fixed by the E-scale fit (see
  llgi_scatfrac_b_factor_target_evaluator). B_scatfrac's sign is
  unconstrained, so both a falling and a rising ScatFrac trend are
  covered; it is restrained toward 0 by
  params.scatfrac_b_factor_restraint_sigma.

  working_selection: flex.bool, True for reflections to include (the
  working set, i.e. NOT r_free_flags -- pass ~r_free_flags.data()).

  The starting ScatFrac_inf is the single-bin ratio-of-sums estimate from
  estimate_llgi_scatfrac; B_scatfrac starts at 0.

  Returns a group_args with .scatfrac (flex.double, one value per input
  reflection, evaluated at ALL reflections regardless of
  working_selection), .target (final fitted LLGI target value on the
  working set, for diagnostics/logging), and .scatfrac_inf/.b_scatfrac
  (floats, for logging/diagnostics).
  """
  if(hybrid is not None):
    # exact wherever possible: a sigmaA-dependent switch would make the
    # fitted objective discontinuous (see llgi_exact.h class hybrid)
    hybrid = hybrid.with_rice_kappa(0.0)
  if(params is None):
    params = llgi_sigmaa_scatfrac_params.extract()
  n_refl = f_eff.size()
  assert working_selection.size() == n_refl
  assert f_calc.size() == n_refl
  assert dobs.size() == n_refl
  assert sigmaa.size() == n_refl
  assert teps.size() == n_refl
  assert resn.size() == n_refl
  assert centric_flags.size() == n_refl
  assert d_star_sq.size() == n_refl
  n_work = working_selection.count(True)
  if(n_work == 0):
    raise RuntimeError(
      "No working-set reflections available for the LLGI ScatFrac fit.")
  initial_scatfrac_inf = flex.mean(estimate_llgi_scatfrac(
    f_calc=f_calc, teps=teps, resn=resn, d_star_sq=d_star_sq,
    centric_flags=centric_flags, scale_factor=scale_factor, n_bins=1))
  # ss = d*^2/4, matching the bulk-solvent k_mask convention
  ss = d_star_sq.as_numpy_array() / 4.0
  evaluator = llgi_scatfrac_b_factor_target_evaluator(
    f_eff=f_eff,
    selection=working_selection,
    f_calc=f_calc,
    dobs=dobs,
    sigmaa=sigmaa,
    teps=teps,
    resn=resn,
    ss=ss,
    centric_flags=centric_flags,
    scale_factor=scale_factor,
    scatfrac_inf_start=initial_scatfrac_inf,
    b_scatfrac_start=0.0,
    restraint_sigma=params.scatfrac_b_factor_restraint_sigma,
    max_iterations=params.max_iterations, hybrid=hybrid)
  scatfrac_inf, b_scatfrac = evaluator.scatfrac_inf_and_b()
  scatfrac = evaluator.scatfrac()
  final_result = ext.llgi_sigmaa_scatfrac_target_and_gradients(
    f_eff=f_eff, selection=working_selection, f_calc=f_calc, dobs=dobs,
    sigmaa=sigmaa, scatfrac=scatfrac, scale_factor=scale_factor,
    teps=teps, resn=resn, centric_flags=centric_flags, hybrid=hybrid)
  return group_args(
    scatfrac=scatfrac, target=final_result.target(),
    scatfrac_inf=scatfrac_inf, b_scatfrac=b_scatfrac)

