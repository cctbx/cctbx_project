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
  sigmaa_curvature_weight = 0.02
    .type = float
    .short_caption = SigmaA spline curvature restraint weight
    .help = "Weight of a light restraint on the second difference of the "\
            "sigmaA B-spline's raw (pre-sigmoid) coefficients, penalising "\
            "curvature of the fitted curve. A sudden collapse of sigmaA "\
            "toward its lower bound at high resolution is not physically "\
            "expected (fit quality should vary smoothly with resolution), "\
            "but can happen where the R-free test set is sparse and the "\
            "LLGI likelihood is nearly flat. This restraint is quadratic "\
            "in the coefficients and vanishes for an exactly log-linear-"\
            "like (constant second difference of z) curve, so it has "\
            "essentially no effect where the data are informative and "\
            "only damps curvature where the likelihood alone cannot "\
            "constrain the fit. See mmtbx.refinement.llgi_sigmaa."\
            "_spline_curvature_penalty_and_gradient. 0 disables it."
  max_iterations = 100
    .type = int
    .expert_level = 3
  scatfrac_curvature_weight = 0.0
    .type = float
    .short_caption = ScatFrac spline curvature restraint weight (EXPERIMENTAL, off by default)
    .help = "Weight of a light restraint on the second difference of the "\
            "ScatFrac B-spline's raw (log-space) coefficients, only used "\
            "when estimate_scatfrac_by_likelihood is True (the empirical "\
            "moment estimator, estimate_llgi_scatfrac, is not restrained "\
            "this way -- it is already a robust binned estimate). Same "\
            "penalty as sigmaa_curvature_weight, restraining z where "\
            "ScatFrac = exp(z) -- here EXACTLY log(ScatFrac), not just "\
            "approximately as for sigmaA's sigmoid, so there is no "\
            "regime-dependence to the restraint's meaning. DEFAULT 0 (OFF), "\
            "UNLIKE sigmaa_curvature_weight: a real 5-macrocycle weight "\
            "sweep on 2G38 (0, 0.001-0.02 coarse, 0.002-0.008 fine) found "\
            "no clean, monotonic safe region -- most weights tested (incl. "\
            "some as small as 0.001) drove at least one macrocycle's own "\
            "ScatFrac-fit objective strongly POSITIVE (worse than the null "\
            "hypothesis, up to +0.48) and gave worse final R-free than no "\
            "restraint at all (e.g. weight=0.02, matching sigmaa_curvature_"\
            "weight's default, gave R-free 0.4250 vs 0.4119 unrestrained); "\
            "only two isolated weights in the whole sweep (0.005, 0.006) "\
            "stayed healthy throughout, with no weights on either side of "\
            "them behaving similarly -- evidence of a genuinely multi-"\
            "modal optimisation landscape (LBFGS landing in a different "\
            "local optimum depending on restraint strength interacting "\
            "with each macrocycle's changing, sometimes marginal, model "\
            "state) rather than a simple over/under-smoothing trade-off. "\
            "Unlike sigmaa_curvature_weight (any weight >= ~0.001 fixed "\
            "the sigmaA collapse cleanly), this restraint is NOT yet safe "\
            "to enable by default; the ScatFrac wobble it was written to "\
            "address (see doc/llgi_target_design.md's ScatFrac restraint "\
            "note) is real but occasional, and leaving this off (as-fit, "\
            "with no curvature restraint at all) is currently the more "\
            "reliable choice. See mmtbx.refinement.llgi_sigmaa."\
            "llgi_scatfrac_target_evaluator and "\
            "_spline_curvature_penalty_and_gradient for the mechanism; "\
            "set > 0 to experiment, at your own risk, until the multi-"\
            "modality above is understood."
  scatfrac_b_factor_restraint_sigma = 10.0
    .type = float
    .short_caption = ScatFrac B-factor restraint sigma (Angstrom^2, b_factor mode only)
    .help = "Weak quadratic restraint on B_scatfrac toward 0, used only "\
            "when scatfrac_model == b_factor. Motivation: sigmaA and "\
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
  scatfrac_model = spline scalar *b_factor
    .type = choice
    .short_caption = ScatFrac(resolution) functional form (Scheme B only)
    .help = "Only used when estimate_scatfrac_by_likelihood is True. "\
            "Chooses how many degrees of freedom the ScatFrac fit gets:\n"\
            "  spline   -- resolution-dependent B-spline curve (see "\
            "llgi_scatfrac_target_evaluator). Historical default for "\
            "Scheme B, but a weight sweep of its curvature restraint "\
            "(scatfrac_curvature_weight) found no safe operating point "\
            "on real data -- most weights made the fit measurably worse "\
            "and occasionally drove it strongly worse than the null "\
            "hypothesis (see doc/llgi_target_design.md's ScatFrac "\
            "restraint note); use with caution.\n"\
            "  scalar   -- ONE ScatFrac value shared by every reflection "\
            "(see llgi_scatfrac_scalar_target_evaluator). No spline "\
            "degrees of freedom for a noisy resolution shell to hijack, "\
            "no curvature to restrain -- robust on real data, matching "\
            "the best hand-tuned spline weight with no tuning. Cannot "\
            "represent genuine resolution dependence in ScatFrac at all.\n"\
            "  b_factor -- DEFAULT. ScatFrac_inf * exp(-B_scatfrac * ss), "\
            "a two-parameter log-linear curve (see llgi_scatfrac_b_factor_"\
            "target_evaluator), directly analogous to an ordinary "\
            "crystallographic B-factor falloff but with B_SCATFRAC'S "\
            "SIGN UNCONSTRAINED (a partial model built from its best-"\
            "ordered components first can show ScatFrac RISING toward "\
            "high resolution, the opposite of the usual falloff). Same "\
            "robustness argument as scalar (a line has no curvature to "\
            "restrain either), while still allowing a genuine monotonic "\
            "resolution trend; the middle ground between spline and "\
            "scalar, and on 2G38 (5 macrocycles) matches the best hand-"\
            "tuned spline weight's R-free with no per-dataset tuning, "\
            "once B_scatfrac's own restraint (scatfrac_b_factor_"\
            "restraint_sigma) is applied -- see doc/llgi_target_design.md's "\
            "ScatFrac restraint note. scatfrac_curvature_weight is "\
            "ignored in both scalar and b_factor modes."
  estimate_scatfrac_by_likelihood = True
    .type = bool
    .short_caption = Fit sigmaA (E-scale) then ScatFrac (F-scale LLGI), instead of a moment estimator
    .help = "DEFAULT True: break the sigmaA/ScatFrac non-identifiability "\
            "the other way round from the historical scheme: fit sigmaA "\
            "FIRST against the E-scale LLGI target (where ScatFrac does "\
            "not appear at all, so there is no degeneracy at that stage "\
            "-- see mmtbx.refinement.llgi_e_bulk_solvent."\
            "estimate_e_sigmaa_fixed_bulk_solvent), then fit ScatFrac as "\
            "a genuine LLGI-likelihood optimisation (F-scale, full "\
            "working set, sigmaA held fixed -- see "\
            "mmtbx.refinement.llgi_sigmaa.estimate_llgi_scatfrac_"\
            "likelihood) instead of the empirical ratio-of-sums moment "\
            "estimator (estimate_llgi_scatfrac). ScatFrac is NOT bounded "\
            "above by 1 in this mode: Feff need not be on absolute scale "\
            "(its scaling assumes 50% solvent content by default), so "\
            "ScatFrac -- a ratio against that scale -- can genuinely "\
            "exceed 1 without indicating a bug. If False, the original "\
            "scheme is used instead: ScatFrac by moment estimator (full "\
            "set), THEN sigmaA by F-scale LLGI (R-free)."
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
        max_iterations=100,
        curvature_weight=0.0):
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
    self.curvature_weight = curvature_weight
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
    penalty, penalty_grad = _spline_curvature_penalty_and_gradient(
      np.array(self.x), self.curvature_weight)
    f += penalty
    g = g + penalty_grad
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
  defaults (only n_sigmaa_coeffs, spline_degree, max_iterations,
  sigmaa_curvature_weight are used here; n_scatfrac_coeffs is used by
  estimate_llgi_scatfrac).

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
    max_iterations=params.max_iterations,
    curvature_weight=params.sigmaa_curvature_weight)
  sigmaa = evaluator.sigmaa()
  final_result = ext.llgi_sigmaa_scatfrac_target_and_gradients(
    f_eff=f_eff, selection=r_free_flags, f_calc=f_calc, dobs=dobs,
    sigmaa=sigmaa, scatfrac=scatfrac, scale_factor=scale_factor,
    teps=teps, resn=resn, centric_flags=centric_flags)
  return group_args(sigmaa=sigmaa, target=final_result.target())

class llgi_scatfrac_target_evaluator(object):
  """ scitbx.lbfgs target evaluator optimising the B-spline coefficients
  of ScatFrac(resolution) against the F-scale LLGI target, summed over
  the FULL (working) reflection set, with sigmaA held fixed -- the
  companion to llgi_sigmaa_target_evaluator, which does the reverse
  (sigmaA free, ScatFrac fixed).

  Design rationale (see doc/llgi_target_design.md's "E-then-F sigmaA/
  ScatFrac estimation" note): sigmaA and ScatFrac are not jointly
  identifiable from the F-scale LLGI target alone (D = Dobs*sigmaA/
  sqrt(ScatFrac) is the only combined quantity the target sees), which is
  why estimate_llgi_scatfrac historically broke the degeneracy with a
  direct empirical (non-likelihood) ScatFrac estimate. This evaluator
  instead breaks the degeneracy the other way: sigmaA is fit FIRST
  against the E-SCALE LLGI target (mmtbx.refinement.llgi_e_bulk_solvent.
  estimate_e_sigmaa/estimate_e_sigmaa_fixed_bulk_solvent), where ScatFrac
  does not appear at all (E-scale Emodel is normalised by sqrt(EPS*
  SigmaP), not by ScatFrac) -- so there is no degeneracy to break at that
  stage. With sigmaA already fixed, ScatFrac can THEN be fit as a genuine
  LLGI-likelihood optimisation (not a moment/ratio estimator) on the full
  working set, using the existing d(target)/d(ScatFrac) derivative
  (ext.llgi_sigmaa_scatfrac_target_and_gradients).

  Reparameterisation: ScatFrac = exp(z), z = design . coeffs -- log-space,
  matching estimate_llgi_scatfrac's own convention, and for the same
  reason: guarantees ScatFrac > 0 automatically with no clipping, and
  (unlike sigmaA's (0.01, 0.99) sigmoid) places NO upper bound on
  ScatFrac. An upper bound of 1 would be wrong here: Feff is not
  guaranteed to be on absolute scale (its scaling assumes 50% solvent
  content by default, which need not hold), so ScatFrac -- defined as a
  ratio against Feff/Resn's absolute scale -- can genuinely exceed 1
  without indicating a bug; see doc/llgi_target_design.md's discussion of
  the observed ScatFrac > 1 values and their likely cause.

  Curvature restraint: the shared second-difference penalty
  (_spline_curvature_penalty_and_gradient, same helper sigmaA's evaluator
  uses) is applied directly to z, exactly as for sigmaA -- but here it is
  restraining log(ScatFrac) EXACTLY, not just approximately in some
  regime: since ScatFrac = exp(z) everywhere (not only near a bound, as
  for sigmaA's sigmoid), z IS log(ScatFrac), so there is no approximation
  to worry about. An earlier version of this evaluator had no restraint
  at all, reasoning that ScatFrac's collapse pathology (the one that
  motivated sigmaA's restraint) had never been observed for ScatFrac and
  that the empirical starting guess (estimate_llgi_scatfrac) should be
  enough on its own -- that reasoning was incomplete: a real 5-macrocycle
  run on 2G38 showed one macrocycle's ScatFrac fit developing a genuine,
  if modest, high-resolution dip not present in neighbouring macrocycles
  (mean-field target sigmaA/ScatFrac wobbling as the model changes cycle
  to cycle, not a bug in the fit itself, but exactly the kind of
  isolated, resolution-localised wobble sigmaA's own restraint was
  designed to damp). See doc/llgi_target_design.md's ScatFrac restraint
  note for the real-data evidence.
  """

  def __init__(self,
        f_eff,
        selection,
        f_calc,
        dobs,
        sigmaa,
        teps,
        resn,
        centric_flags,
        scale_factor,
        scatfrac_design,
        n_scatfrac_coeffs,
        scatfrac_start,
        max_iterations=100,
        curvature_weight=0.0):
    self.f_eff = f_eff
    self.selection = selection
    self.f_calc = f_calc
    self.dobs = dobs
    self.sigmaa = sigmaa
    self.teps = teps
    self.resn = resn
    self.centric_flags = centric_flags
    self.scale_factor = scale_factor
    self.scatfrac_design = scatfrac_design  # numpy array, (n_refl, n_coeffs)
    self.n_scatfrac_coeffs = n_scatfrac_coeffs
    self.curvature_weight = curvature_weight
    # Start from log(scatfrac_start) (e.g. estimate_llgi_scatfrac's own
    # empirical estimate), least-squares-fit onto the spline basis at the
    # SAME per-reflection points this evaluator will evaluate at -- not a
    # neutral z=0 start (unlike sigmaA): ScatFrac already has a good
    # starting guess available, and starting LBFGS from it is both faster
    # and less likely to wander into a degenerate/implausible region than
    # starting flat.
    import numpy as np
    log_start = np.log(np.asarray(scatfrac_start, dtype=float))
    coeffs0, _res, _rank, _sv = np.linalg.lstsq(
      scatfrac_design, log_start, rcond=None)
    self.x = flex.double(coeffs0)
    term_parameters = scitbx.lbfgs.termination_parameters(
      max_iterations=max_iterations)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=term_parameters)

  def _current_scatfrac(self):
    import numpy as np
    coeffs = np.array(self.x)
    z = self.scatfrac_design.dot(coeffs)
    scatfrac = np.exp(z)
    return scatfrac, scatfrac  # d(exp(z))/dz == exp(z) == scatfrac itself

  def compute_functional_and_gradients(self):
    import numpy as np
    scatfrac, dscatfrac_dz = self._current_scatfrac()
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
      centric_flags=self.centric_flags)
    f = result.target()
    d_target_by_dscatfrac = np.array(result.d_target_by_dscatfrac())
    # Chain rule: d(target)/d(coeffs) = design.T @ (d_target * dscatfrac_dz),
    # since z_i = design[i,:] . coeffs is linear in coeffs.
    g = self.scatfrac_design.T.dot(d_target_by_dscatfrac * dscatfrac_dz)
    penalty, penalty_grad = _spline_curvature_penalty_and_gradient(
      np.array(self.x), self.curvature_weight)
    f += penalty
    g = g + penalty_grad
    return f, flex.double(g)

  def scatfrac(self):
    """ Return the fitted ScatFrac array at the current (final, once the
    minimizer has run) coefficient values. """
    scatfrac, _ = self._current_scatfrac()
    return flex.double(scatfrac)

class llgi_scatfrac_scalar_target_evaluator(object):
  """ scitbx.lbfgs target evaluator optimising a SINGLE ScatFrac value
  (not a resolution-dependent spline) against the F-scale LLGI target,
  summed over the full working set, sigmaA held fixed -- the scalar
  counterpart of llgi_scatfrac_target_evaluator.

  Motivation (see doc/llgi_target_design.md's ScatFrac restraint note):
  a weight sweep of the spline curvature restraint on real 2G38 data
  found no safe operating point -- most weights tested made the fit
  measurably worse and occasionally drove the per-macrocycle objective
  strongly positive, in a non-monotonic, multi-modal pattern (isolated
  safe weights surrounded by unsafe ones on both sides), evidence that
  the spline's extra degrees of freedom let LBFGS chase noise in sparse
  resolution shells into a genuinely different, worse local optimum from
  one macrocycle to the next. A single scalar has no such freedom: no
  per-region basis functions for a noisy shell to hijack, no curvature to
  restrain (a constant has none by construction), and a much smaller,
  better-conditioned optimisation (one parameter, not
  n_scatfrac_coeffs). This is a coarser model of ScatFrac(resolution) --
  it cannot represent genuine resolution dependence in the fraction of
  scattering accounted for -- but is offered as a more robust default
  than the spline for exactly the reason the spline became unsafe.

  Same log-space reparameterisation as the spline version (ScatFrac =
  exp(z), NO upper bound of 1 -- see llgi_scatfrac_target_evaluator's
  docstring for why an upper bound would be wrong), started from the
  weighted mean of estimate_llgi_scatfrac's own empirical per-bin
  estimate rather than a neutral guess, for the same reason the spline
  version starts from that estimator's output.
  """

  def __init__(self,
        f_eff,
        selection,
        f_calc,
        dobs,
        sigmaa,
        teps,
        resn,
        centric_flags,
        scale_factor,
        scatfrac_start,
        max_iterations=100):
    self.f_eff = f_eff
    self.selection = selection
    self.f_calc = f_calc
    self.dobs = dobs
    self.sigmaa = sigmaa
    self.teps = teps
    self.resn = resn
    self.centric_flags = centric_flags
    self.scale_factor = scale_factor
    self.n_refl = f_eff.size()
    import numpy as np
    import math
    z0 = math.log(float(scatfrac_start))
    self.x = flex.double([z0])
    term_parameters = scitbx.lbfgs.termination_parameters(
      max_iterations=max_iterations)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=term_parameters)

  def _current_scatfrac(self):
    import math
    z = self.x[0]
    scatfrac = math.exp(z)
    return scatfrac, scatfrac  # d(exp(z))/dz == exp(z) == scatfrac itself

  def compute_functional_and_gradients(self):
    import numpy as np
    scatfrac, dscatfrac_dz = self._current_scatfrac()
    result = ext.llgi_sigmaa_scatfrac_target_and_gradients(
      f_eff=self.f_eff,
      selection=self.selection,
      f_calc=self.f_calc,
      dobs=self.dobs,
      sigmaa=self.sigmaa,
      scatfrac=flex.double(self.n_refl, scatfrac),
      scale_factor=self.scale_factor,
      teps=self.teps,
      resn=self.resn,
      centric_flags=self.centric_flags)
    f = result.target()
    d_target_by_dscatfrac = np.array(result.d_target_by_dscatfrac())
    # Chain rule: every reflection's ScatFrac is the SAME scalar, so
    # d(target)/dz = sum_i d_target_by_dscatfrac[i] * dscatfrac_dz
    # (unlike the spline case, no design matrix -- every reflection
    # contributes to the one shared parameter directly).
    g = float(np.sum(d_target_by_dscatfrac)) * dscatfrac_dz
    return f, flex.double([g])

  def scatfrac(self):
    """ Return the fitted scalar ScatFrac value, broadcast to one entry
    per input reflection (matching llgi_scatfrac_target_evaluator's
    return shape, so callers do not need to special-case this mode). """
    scatfrac, _ = self._current_scatfrac()
    return flex.double(self.n_refl, scatfrac)

class llgi_scatfrac_b_factor_target_evaluator(object):
  """ scitbx.lbfgs target evaluator optimising a TWO-PARAMETER ScatFrac
  curve, ScatFrac_inf * exp(-B_scatfrac * ss) (ss = d*^2/4, the same
  sin(theta)/lambda squared convention as the bulk-solvent k_mask formula
  -- see mmtbx.refinement.llgi_e_bulk_solvent.k_mask_and_gradients -- so
  B_scatfrac is directly comparable to an ordinary crystallographic
  B-factor, positive B meaning ScatFrac falls off toward high resolution
  as usual), against the F-scale LLGI target, sigmaA held fixed -- a
  middle ground between llgi_scatfrac_target_evaluator (the unrestrained
  spline, unsafe -- see its docstring) and llgi_scatfrac_scalar_target_
  evaluator (one value, no resolution dependence at all).

  Motivation: the scalar mode's own docstring notes it cannot represent
  genuine resolution dependence in the fraction of scattering accounted
  for. That dependence is physically plausible in more than one
  direction: a partial model built from its best-ordered components
  first (as coordinate refinement of an initially poor model often
  proceeds) would have ScatFrac RISE toward high resolution (the model
  accounts for a larger share of the weak, well-ordered high-resolution
  data than of the disordered low-resolution features it has not yet
  captured) -- the opposite sign from the more familiar case of ScatFrac
  falling with resolution. A single positive-only B-factor-style
  parameter B_scatfrac (sign unconstrained here, unlike an ordinary ADP)
  covers both directions with only one more parameter than the scalar
  fit, still far short of the spline's degrees of freedom.

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

  Same log-space reparameterisation as the other two ScatFrac modes
  (ScatFrac = exp(z), NO upper bound of 1). No curvature restraint is
  relevant here either: a straight line in log-ScatFrac-vs-ss space has
  no curvature to restrain by construction (same reasoning as the scalar
  mode), so this mode carries the same robustness argument that motivated
  the scalar mode over the spline.

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
  through the F-scale LLGI target (see llgi_sigmaa_scatfrac_params.
  estimate_scatfrac_by_likelihood's docstring), so a B_scatfrac far
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
        max_iterations=100):
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
      centric_flags=self.centric_flags)
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
      params=None):
  """ Fit ScatFrac(resolution) as a B-spline curve AGAINST THE F-SCALE
  LLGI TARGET (not the moment/ratio estimator in estimate_llgi_scatfrac),
  on the full working set (all reflections minus R-free -- see the design
  note's "working set, not all data" rationale, same as the E-scale bulk-
  solvent fit's own data split), with sigmaA already fixed (see
  llgi_scatfrac_target_evaluator's docstring for why fixing sigmaA FIRST,
  via the E-scale target where ScatFrac does not appear, breaks the
  sigmaA/ScatFrac degeneracy without needing a non-likelihood estimator).

  working_selection: flex.bool, True for reflections to include (the
  working set, i.e. NOT r_free_flags -- pass ~r_free_flags.data()).

  Uses estimate_llgi_scatfrac's own empirical ratio-of-sums estimate only
  as the LBFGS starting point (see llgi_scatfrac_target_evaluator), not
  as the returned answer. A light curvature restraint (params.scatfrac_
  curvature_weight, same second-difference penalty as sigmaA's) damps
  resolution-localised wobbles in the fitted curve -- see llgi_scatfrac_
  target_evaluator's docstring for the real-data motivation, and its
  DEFAULT-OFF status (a weight sweep on real data found no safe operating
  point -- see llgi_scatfrac_scalar_target_evaluator's docstring for the
  more robust alternative this function can dispatch to).

  params.scatfrac_model selects the functional form (see that phil
  parameter's help for the full rationale): "spline" (historical
  default, unsafe -- see llgi_scatfrac_target_evaluator), "scalar" (one
  value, see llgi_scatfrac_scalar_target_evaluator), or "b_factor" (a
  two-parameter ScatFrac_inf*exp(-B_scatfrac*ss) curve, see llgi_
  scatfrac_b_factor_target_evaluator -- B_scatfrac's sign is
  unconstrained, so this covers both a falling AND a rising ScatFrac
  trend, restrained toward 0 by params.scatfrac_b_factor_restraint_sigma
  -- see llgi_scatfrac_b_factor_target_evaluator's docstring). params.
  scatfrac_curvature_weight is ignored except in "spline" mode (neither
  "scalar" nor "b_factor" has curvature to restrain).

  Returns a group_args with .scatfrac (flex.double, one value per input
  reflection, evaluated at ALL reflections regardless of
  working_selection) and .target (final fitted LLGI target value on the
  working set, for diagnostics/logging). In "b_factor" mode, also
  .scatfrac_inf and .b_scatfrac (floats, for logging/diagnostics).
  """
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
  scatfrac_inf = None
  b_scatfrac = None
  if(params.scatfrac_model == "scalar"):
    # A single-bin call to estimate_llgi_scatfrac still returns a
    # per-reflection flex.double (the spline evaluated everywhere, from
    # fitting to ONE bin's point estimate) rather than a bare scalar --
    # take its mean as the scalar starting point (harmless even if the
    # single-bin fit is not perfectly constant across reflections, since
    # this is only an LBFGS starting guess, not the returned answer).
    initial_scatfrac_scalar = flex.mean(estimate_llgi_scatfrac(
      f_calc=f_calc, teps=teps, resn=resn, d_star_sq=d_star_sq,
      centric_flags=centric_flags, scale_factor=scale_factor,
      n_coeffs=params.n_scatfrac_coeffs,
      spline_degree=params.spline_degree, n_bins=1))
    evaluator = llgi_scatfrac_scalar_target_evaluator(
      f_eff=f_eff,
      selection=working_selection,
      f_calc=f_calc,
      dobs=dobs,
      sigmaa=sigmaa,
      teps=teps,
      resn=resn,
      centric_flags=centric_flags,
      scale_factor=scale_factor,
      scatfrac_start=initial_scatfrac_scalar,
      max_iterations=params.max_iterations)
  elif(params.scatfrac_model == "b_factor"):
    # Same single-bin empirical estimate as the scalar mode's own
    # starting point; B_scatfrac starts at 0 (flat), a neutral guess
    # that does not presuppose a falling or rising trend.
    initial_scatfrac_inf = flex.mean(estimate_llgi_scatfrac(
      f_calc=f_calc, teps=teps, resn=resn, d_star_sq=d_star_sq,
      centric_flags=centric_flags, scale_factor=scale_factor,
      n_coeffs=params.n_scatfrac_coeffs,
      spline_degree=params.spline_degree, n_bins=1))
    # ss = d*^2/4 (sin(theta)/lambda squared), matching the bulk-solvent
    # k_mask formula's own convention (see mmtbx.refinement.
    # llgi_e_bulk_solvent.k_mask_and_gradients/ss_from_f_obs) -- computed
    # directly here rather than imported, to avoid a cross-module
    # dependency for one line.
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
      max_iterations=params.max_iterations)
    scatfrac_inf, b_scatfrac = evaluator.scatfrac_inf_and_b()
  else:
    assert params.scatfrac_model == "spline", params.scatfrac_model
    initial_scatfrac = estimate_llgi_scatfrac(
      f_calc=f_calc, teps=teps, resn=resn, d_star_sq=d_star_sq,
      centric_flags=centric_flags, scale_factor=scale_factor,
      n_coeffs=params.n_scatfrac_coeffs, spline_degree=params.spline_degree,
      n_bins=params.n_scatfrac_bins)
    scatfrac_design = _b_spline_design_matrix(
      d_star_sq.as_numpy_array(), params.n_scatfrac_coeffs,
      params.spline_degree)
    evaluator = llgi_scatfrac_target_evaluator(
      f_eff=f_eff,
      selection=working_selection,
      f_calc=f_calc,
      dobs=dobs,
      sigmaa=sigmaa,
      teps=teps,
      resn=resn,
      centric_flags=centric_flags,
      scale_factor=scale_factor,
      scatfrac_design=scatfrac_design,
      n_scatfrac_coeffs=params.n_scatfrac_coeffs,
      scatfrac_start=initial_scatfrac,
      max_iterations=params.max_iterations,
      curvature_weight=params.scatfrac_curvature_weight)
  scatfrac = evaluator.scatfrac()
  final_result = ext.llgi_sigmaa_scatfrac_target_and_gradients(
    f_eff=f_eff, selection=working_selection, f_calc=f_calc, dobs=dobs,
    sigmaa=sigmaa, scatfrac=scatfrac, scale_factor=scale_factor,
    teps=teps, resn=resn, centric_flags=centric_flags)
  return group_args(
    scatfrac=scatfrac, target=final_result.target(),
    scatfrac_inf=scatfrac_inf, b_scatfrac=b_scatfrac)

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
