from __future__ import absolute_import, division, print_function
import numpy as np
import scitbx.minimizers
from cctbx.array_family import flex
from libtbx import group_args
import iotbx.phil
import mmtbx.refinement.llgi_e_dmodel as dmodel
import mmtbx.refinement.llgi_e_dmodel_target as target

""" L-BFGS-B fit of the physically-motivated D_model(s; theta)
parametrization against the E-scale LLGI likelihood (see
doc/llgi_target_design.md sec. 6.4), a drop-in replacement for
mmtbx.refinement.llgi_e_bulk_solvent.estimate_e_sigmaa's
B-spline-over-sigmoid fit in the Stage-1/Stage-2 inner loop
(run_inner_loop), with the same restriction to the R-free/test set and
the same Emodel-held-fixed convention.

Builds on llgi_e_dmodel.py (D_model value/gradient/Hessian),
llgi_e_likelihood.py (per-reflection likelihood) and
llgi_e_dmodel_target.py (chain-rule combination into theta-space
LL/gradient/Hessian).
"""

llgi_e_dmodel_params = iotbx.phil.parse("""\
  n_gaussian_terms = 2
    .type = int
    .short_caption = Number of coordinate-error Gaussian terms (K)
    .help = "Number of positive Gaussian-decay terms in D_model(s; "\
            "theta) (sigmaA_model_handoff.md sec. 2), each standing in "\
            "for a discretized component of the distribution of "\
            "coordinate-error B-factors across the model. K=2 or 3 is "\
            "expected to be sufficient (design doc sec. 6.4)."
  max_iterations = 200
    .type = int
    .expert_level = 3
  include_constant_term = False
    .type = bool
    .short_caption = Add a constant (B=0) term to the D_model ladder
    .help = "Prepend a B=0 rung (a resolution-independent amplitude) to "\
            "the n_gaussian_terms log-spaced B_k ladder. The ladder's "\
            "lowest decaying rung (4*d_min^2) has already fallen to 1/e "\
            "at d_min, so without this term D_model cannot stay flat "\
            "at high resolution, as sigmaA does for a well-refined "\
            "model (seen on 9RRL). With it, D_model need not fall to 0 "\
            "at infinite resolution."
  a_k_smoothness_weight = 0.1
    .type = float
    .short_caption = Smoothness penalty weight across the fixed B_k ladder
    .help = "Weight of a roughness (second-difference) penalty across "\
            "the K coordinate-error amplitudes a_k, in the ORDER of "\
            "the fixed b_k_grid ladder they multiply (llgi_e_dmodel_"\
            "fit.default_b_k_grid) -- penalty = weight * sum_k "\
            "(a_{k+1} - 2*a_k + a_{k-1})^2 over interior rungs (a no-op "\
            "for K<3, where no interior rung exists). Exists for the "\
            "same reason a_k_smoothness_penalty_and_gradient's own "\
            "docstring gives in full: theta no longer includes B_k at "\
            "all (design doc sec. 6.4's well-posedness addendum -- B_k "\
            "is now a FIXED ladder, not fitted, precisely because the "\
            "earlier free-B_k parametrization was an ill-posed, nearly-"\
            "collinear nonlinear fit that a real synthetic-data test "\
            "(tst_llgi_e_dmodel_fit.py) drove to a degenerate a_k->0/"\
            "B_k->infinity solution once D_model's high-resolution "\
            "asymptote was corrected to genuinely reach 0), so the "\
            "coordinate-error-varies-smoothly-with-B-factor physical "\
            "picture this restraint encodes is now expressed directly "\
            "as smoothness across NEIGHBOURING a_k on that fixed "\
            "ladder, rather than as a restraint on B_k itself (no "\
            "longer a free parameter to restrain). 0 (or None) "\
            "disables it -- becomes an unrestrained non-negative-least-"\
            "squares-like fit of a_k against the ladder, still well-"\
            "posed (unlike the old free-B_k form) but no longer "\
            "favouring a smooth coefficient profile over a jagged one."
""")

def default_b_k_grid(k, s2):
  """ Fixed, log-spaced ladder of K coordinate-error decay constants
  B_1..B_K (design doc sec. 6.4's well-posedness addendum -- see
  llgi_e_dmodel.py's own module docstring for why B_k is a fixed grid
  rather than a fitted parameter). Spans the resolution range actually
  present in the fitted (test-set) reflections, s2 = d_star_sq/4: from
  roughly 1/s2_max (a term that has already decayed to ~1/e by the
  data's own highest-resolution reflection -- any B_k smaller than this
  is indistinguishable from a constant over the whole fitted range) to
  roughly 1/s2_min (a term that has already decayed to ~1/e by the
  data's own LOWEST-resolution reflection -- any B_k larger than this
  contributes essentially nothing anywhere in the fitted range), evenly
  spaced in log(B_k) (matching how coordinate-error B-factors are
  naturally compared, and how sigmaA_model_handoff.md sec. 2 frames the
  K terms as "a discretized component of the distribution").

  k: number of ladder rungs (K>=1). s2: the fitted reflections' own s^2
  = d_star_sq/4 values (array). Falls back to a fixed, physically
  plausible range (10 to 300, design doc sec. 2's own realistic B_k
  span) if s2 is degenerate (empty or a single repeated value).

  Returns a 1D numpy array of length k, ascending.
  """
  if(k <= 0):
    return np.zeros(0, dtype=float)
  s2 = np.asarray(s2, dtype=float)
  s2_pos = s2[s2 > 0.0]
  if(s2_pos.size == 0 or float(np.min(s2_pos)) == float(np.max(s2_pos))):
    b_lo, b_hi = 10.0, 300.0
  else:
    s2_min, s2_max = float(np.min(s2_pos)), float(np.max(s2_pos))
    b_lo, b_hi = 1.0 / s2_max, 1.0 / s2_min
  if(k == 1):
    return np.array([np.sqrt(b_lo * b_hi)], dtype=float)
  return np.exp(np.linspace(np.log(b_lo), np.log(b_hi), k))

def a_k_smoothness_penalty_and_gradient(theta, weight, first=0):
  """ Roughness penalty discouraging a jagged a_k profile across the
  fixed b_k_grid ladder (llgi_e_dmodel_params.a_k_smoothness_weight's
  own docstring has the full motivation: this replaces the earlier
  free-B_k restraint now that B_k is fixed, not fitted). Standard
  second-difference (discrete-Laplacian) Tikhonov form:

    penalty = weight * sum_{k=1}^{K-2} (a_{k+1} - 2*a_k + a_{k-1})^2

  (0-indexed interior rungs only -- a no-op for K<3, where there is no
  interior rung to take a second difference at). This is a genuine
  smoothness prior on the FITTED coefficients themselves (not a proxy
  restraint on an already-fixed quantity), directly matching the
  physical picture of coordinate error varying smoothly across a
  continuum of atomic B-factors (adjacent ladder rungs are adjacent
  points on that continuum).

  theta: natural-space parameter vector. weight: penalty weight
  (a_k_smoothness_weight); <= 0 or None returns a no-op. first: number of
  leading rungs left out of the penalty -- the B=0 constant rung added
  by include_constant_term is not a point on the coordinate-error
  continuum (log B = -infinity), so smoothness across it would tie the
  constant amplitude to the decaying ones (a_0 ~ 2*a_1 - a_2), letting
  the nested model fit worse than plain D_model.

  Returns (penalty, d(penalty)/d(theta)), penalty a plain float, the
  gradient a numpy array the same shape as theta (nonzero only in the
  a_k slots).
  """
  theta = np.asarray(theta, dtype=float)
  grad = np.zeros_like(theta)
  if(weight is None or weight <= 0):
    return 0.0, grad
  a, _, _ = dmodel.unpack_theta(theta)
  a = a[first:]
  k = a.size
  if(k < 3):
    return 0.0, grad
  d2 = a[2:] - 2.0 * a[1:-1] + a[:-2]  # shape (k-2,), d2[i] = 2nd diff at rung i+1
  penalty = weight * float(np.sum(d2 * d2))
  # d(penalty)/da_j = weight * sum_i 2*d2[i] * d(d2[i])/da_j; each d2[i]
  # touches rungs i, i+1, i+2 with coefficients +1, -2, +1 respectively.
  g = np.zeros(k, dtype=float)
  g[:-2] += 2.0 * weight * d2         # d2[i]'s +1 coefficient on rung i
  g[1:-1] += -4.0 * weight * d2       # d2[i]'s -2 coefficient on rung i+1
  g[2:] += 2.0 * weight * d2          # d2[i]'s +1 coefficient on rung i+2
  grad[first:first + k] = g
  return penalty, grad

class d_model_target_evaluator(object):
  """ L-BFGS-B fit of D_model(s; theta) against the E-scale LLGI target,
  summed over the R-free/test set only (same restriction as
  mmtbx.refinement.llgi_e_bulk_solvent.e_sigmaa_target_evaluator), with
  Emodel (hence the bulk-solvent model) held fixed.

  B_defect is fixed at b_sol_anchor (the bulk-solvent fit's B_sol) when
  one is given, and only a_1..a_K and b are fitted. Leaving it free
  (even under a restraint) lets the fit trade the defect term against a
  coordinate-error term with a similar decay -- solutions such as
  a_1 = 233, b = 232 -- and makes the result depend on the starting
  point, for very little gain in likelihood. Without an anchor B_defect
  is fitted too.

  The fit works in natural coordinates with bounds (0 <= a_k, b <=
  amplitude_max; b_defect_min <= B_defect <= b_defect_max). An earlier version optimised q = ln(theta)
  with unbounded L-BFGS, which makes 0 an absorbing boundary (d/dq ->
  0 as theta -> 0): once b or an a_k headed towards 0 it could not come
  back, and the fit converged to whichever such corner it fell into
  first (on 2G38 after 5 cycles, LLGI 156.7 instead of 180.8 with b
  stuck at 0).
  """

  b_defect_min = 0.1
  b_defect_max = 1.e4
  # Upper bound for a_k and b: tanh(D_raw) is saturated (> 0.9999) well
  # before D_raw reaches this, so larger values only drift along flat
  # directions (e.g. the largest-B rung, which matters only at low
  # resolution, where the curve is already saturated).
  amplitude_max = 100.

  def __init__(self,
        e_eff, r_free_flags, e_model, dobs, centric_flags, d_star_sq,
        n_gaussian_terms=2, theta_start=None, max_iterations=200,
        b_sol_anchor=None,
        a_k_smoothness_weight=0.1, b_k_grid=None,
        include_constant_term=False,
      hybrid=None):
    n_refl = e_eff.size()
    assert r_free_flags.size() == n_refl
    assert e_model.size() == n_refl
    assert dobs.size() == n_refl
    assert centric_flags.size() == n_refl
    assert d_star_sq.size() == n_refl
    test_sel = np.array(r_free_flags, dtype=bool)
    if(not np.any(test_sel)):
      raise RuntimeError(
        "d_model_target_evaluator: no R-free/test-set reflections "
        "available for the D_model(s) LLGI sigmaA fit.")
    self.s2 = np.array(d_star_sq, dtype=float)[test_sel] / 4.0
    self.e_eff = np.array(e_eff, dtype=float)[test_sel]
    self.e_c = np.array(e_model, dtype=float)[test_sel]
    self.dobs = np.array(dobs, dtype=float)[test_sel]
    self.centric_flags = np.array(centric_flags, dtype=bool)[test_sel]
    self.hybrid = None
    if(hybrid is not None):
      self.hybrid = hybrid.select(flex.bool(test_sel.tolist()))
    # LL and gradient are per-reflection means, matching ext.
    # llgi_e_sigmaa_target_and_gradients (target() is divided by
    # n_selected), so .final_target is comparable with the spline path's.
    self.n_test = int(np.sum(test_sel))
    self.n_gaussian_terms = n_gaussian_terms
    self.a_k_smoothness_weight = a_k_smoothness_weight
    if(b_k_grid is None):
      b_k_grid = default_b_k_grid(n_gaussian_terms, self.s2)
      if(include_constant_term):
        b_k_grid = np.concatenate([[0.0], b_k_grid])
    self.b_k_grid = np.asarray(b_k_grid, dtype=float)
    n_terms = n_gaussian_terms + int(include_constant_term)
    assert self.b_k_grid.shape == (n_terms,)
    # Leading B=0 (constant) rungs are excluded from the smoothness penalty
    self.n_constant_rungs = int(np.sum(self.b_k_grid <= 0))
    assert np.all(self.b_k_grid[:self.n_constant_rungs] <= 0)
    self.b_defect_fixed = None
    if(b_sol_anchor is not None and b_sol_anchor > 0):
      self.b_defect_fixed = float(b_sol_anchor)
    self.final_target = None

    if(theta_start is None):
      theta_start = self._default_theta_start(n_terms)
    theta_start = np.asarray(theta_start, dtype=float)
    assert theta_start.size == n_terms + 2
    # Fitted parameters: all of theta, or all but B_defect when it is fixed
    self.n_fit = n_terms + 2 - int(self.b_defect_fixed is not None)
    lower = np.zeros(self.n_fit)
    upper = np.full(self.n_fit, self.amplitude_max)
    if(self.b_defect_fixed is None):
      lower[-1] = self.b_defect_min
      upper[-1] = self.b_defect_max
    self.x = flex.double(np.clip(theta_start[:self.n_fit], lower, upper))
    self.bound_flags = flex.int(self.n_fit, 2)  # lower and upper bounds
    self.lower_bound = flex.double(lower)
    self.upper_bound = flex.double(upper)
    self.minimizer = scitbx.minimizers.lbfgs(
      mode="lbfgsb", calculator=self, max_iterations=max_iterations)
    self.update(self.minimizer.x)

  @staticmethod
  def _default_theta_start(k):
    """ Neutral starting theta: each a_k = 0.5/K, a small defect term
    (b = 0.05) and B_defect = 40 (used only when B_defect is fitted).
    """
    theta = np.empty(k + 2, dtype=float)
    if(k > 0):
      theta[:k] = 0.5 / k
    theta[-2] = 0.05
    theta[-1] = 40.0
    return theta

  def _full_theta(self, theta_fit):
    theta_fit = np.asarray(theta_fit, dtype=float)
    if(self.b_defect_fixed is None): return theta_fit
    return np.concatenate([theta_fit, [self.b_defect_fixed]])

  # calculator interface for scitbx.minimizers.lbfgs
  def update(self, x):
    self.x = x
    theta = self._full_theta(np.array(x))
    ll, grad_p, _ = target.total_ll_gradient_hessian(
      theta, self.s2, self.e_eff, self.e_c, self.dobs,
      self.centric_flags, self.b_k_grid, hybrid=self.hybrid)
    smooth_penalty, smooth_grad_p = a_k_smoothness_penalty_and_gradient(
      theta, self.a_k_smoothness_weight, first=self.n_constant_rungs)
    # Minimize-me convention (as llgi_e.h's target_one_h): f = -LL/n plus
    # the smoothness penalty, which is already in minimize-me form.
    n = self.n_fit
    self._f = -ll / self.n_test + smooth_penalty
    self._g = -grad_p[:n] / self.n_test + smooth_grad_p[:n]
    self.final_target = self._f

  def target(self):
    return self._f

  def gradients(self):
    return flex.double(self._g)

  def theta(self):
    return self._full_theta(np.array(self.x))

def estimate_d_model_sigmaa(e_eff, r_free_flags, e_model, dobs,
      centric_flags, d_star_sq, n_gaussian_terms=2, max_iterations=200,
      b_sol_anchor=None,
      a_k_smoothness_weight=0.1, theta_start=None, b_k_grid=None,
      include_constant_term=False, hybrid=None):
  """ Fit D_model(s; theta) against the E-scale LLGI target, restricted
  to the R-free/test set, Emodel held fixed -- drop-in replacement for
  mmtbx.refinement.llgi_e_bulk_solvent.estimate_e_sigmaa in the Stage-1/
  Stage-2 inner loop (run_inner_loop), same call-site role. Evaluates
  the fitted curve at every reflection (working set included, unlike
  the fit itself, exactly mirroring estimate_e_sigmaa's own contract).

  b_k_grid: fixed ladder of B_1..B_K decay constants (length
  n_gaussian_terms, plus one if include_constant_term); None (the
  default) derives it from the fitted (test-set) reflections' own
  resolution range via default_b_k_grid, with a leading B=0 rung if
  include_constant_term. Exposed as its own argument mainly for tests/
  diagnostics that need a reproducible, data-independent ladder --
  ordinary callers should leave it None.

  Returns a group_args with .sigmaa (flex.double, D_model(s_h; theta)
  evaluated at every input reflection, one value per input reflection
  -- the "sigmaa" name kept for drop-in compatibility with
  estimate_e_sigmaa's own return contract, even though the underlying
  parametrization is now the physically-motivated D_model, not a
  spline), .theta (flex.double, the converged natural-space parameter
  vector -- [a_1..a_K, b, B_defect], NOT including B_1..B_K, which are
  fixed and available as .b_k_grid instead -- for logging/diagnostics/
  passing to a subsequent macrocycle as theta_start), .b_k_grid (the
  fixed B_k ladder actually used), .target (final fitted LLGI target
  value on the test set, minimize-me convention, for diagnostics/
  logging matching estimate_e_sigmaa's own .target).
  """
  n_refl = e_eff.size()
  evaluator = d_model_target_evaluator(
    e_eff=e_eff, r_free_flags=r_free_flags, e_model=e_model, dobs=dobs,
    centric_flags=centric_flags, d_star_sq=d_star_sq,
    n_gaussian_terms=n_gaussian_terms, theta_start=theta_start,
    max_iterations=max_iterations,
    b_sol_anchor=b_sol_anchor,
    a_k_smoothness_weight=a_k_smoothness_weight,
    b_k_grid=b_k_grid,
    include_constant_term=include_constant_term,
    hybrid=hybrid)
  theta = evaluator.theta()
  b_k_grid_used = evaluator.b_k_grid
  s2_all = np.array(d_star_sq, dtype=float) / 4.0
  sigmaa_all = dmodel.d_model(s2_all, theta, b_k_grid_used)
  d_star_sq_np = np.asarray(d_star_sq, dtype=float)
  x_range = (float(d_star_sq_np.min()), float(d_star_sq_np.max()))

  def evaluate_at(d_star_sq_new, x_range=x_range):
    """ Evaluate this SAME fitted D_model(s; theta) curve at arbitrary
    new d_star_sq values -- e.g. mmtbx.map_tools.
    model_missing_reflections_llgi's missing-reflection fill (same
    role as e_sigmaa_target_evaluator.evaluate_at, kept for drop-in
    compatibility with that caller's contract). Unlike the B-spline
    fit, D_model(s; theta) is a closed-form sum of Gaussians, well-
    defined at ANY s -- no range-clamping is needed (x_range is
    accepted, for interface compatibility only, and otherwise
    ignored), so this simply evaluates the fitted curve directly, even
    for a missing reflection at a resolution outside the range the fit
    itself was built against.
    """
    s2_new = np.asarray(d_star_sq_new, dtype=float) / 4.0
    return flex.double(dmodel.d_model(s2_new, theta, b_k_grid_used))

  return group_args(
    sigmaa=flex.double(sigmaa_all),
    theta=flex.double(theta),
    b_k_grid=flex.double(b_k_grid_used),
    target=evaluator.final_target,
    x_range=x_range,
    evaluate_at=evaluate_at)
