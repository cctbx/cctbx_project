from __future__ import absolute_import, division, print_function
import numpy as np
import scitbx.lbfgs
from cctbx.array_family import flex
from libtbx import group_args
import iotbx.phil
import mmtbx.refinement.llgi_e_dmodel as dmodel
import mmtbx.refinement.llgi_e_dmodel_target as target
import mmtbx.refinement.llgi_e_dmodel_reparam as reparam

""" L-BFGS fit of the physically-motivated D_model(s; theta)
parametrization against the E-scale LLGI likelihood (Section 7,
implementation step 4 -- see doc/llgi_target_design.md sec. 6.4), a
drop-in replacement for mmtbx.refinement.llgi_e_bulk_solvent.
estimate_e_sigmaa's B-spline-over-sigmoid fit in the Stage-1/Stage-2
inner loop (run_inner_loop), same restriction to the R-free/test set,
same Emodel-held-fixed convention.

Builds on llgi_e_dmodel.py (D_model value/gradient/Hessian),
llgi_e_likelihood.py (corrected symmetric-form per-reflection
likelihood), llgi_e_dmodel_target.py (chain-rule combination into
theta-space LL/gradient/Hessian), and llgi_e_dmodel_reparam.py (log/exp
unconstrained reparametrization, INCLUDING the corrected Hessian chain
rule) -- all pure numpy, independently finite-difference-verified. This
module is the first place scitbx.lbfgs/flex machinery and the actual
sign-flip convention (minimize-me, matching every other evaluator class
in this file's sibling module llgi_e_bulk_solvent.py) enter the
picture.
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
  b_sol_restraint_sigma = 0.6931471805599453
    .type = float
    .short_caption = B_defect restraint width in ln(B), default ln(2)
    .help = "Width (in ln(B), i.e. natural-log space) of a one-"\
            "directional soft restraint pulling B_defect toward the "\
            "bulk-solvent fit's own B_sol (log-linear point estimate "\
            "of the per-bin k_mask curve -- mmtbx.refinement."\
            "llgi_e_bulk_solvent._log_linear_k_sol_b_sol), NOT the "\
            "reverse: B_sol is treated as a fixed anchor from the "\
            "(better-conditioned, already-converged) bulk-solvent "\
            "step, B_defect is pulled toward it, not vice versa (see "\
            "doc/llgi_target_design.md sec. 6.4). Default ln(2) "\
            "restrains B_defect and B_sol to within roughly a factor "\
            "of 2 of each other at 1 sigma -- a deliberately weak "\
            "starting value, meant to be tested against 0/None "\
            "(unrestrained) to check whether the two should really be "\
            "expected to agree this closely. B_sol and B_defect are "\
            "physically related but NOT identical quantities (B_sol: "\
            "resolution-decay of the bulk-solvent MEAN amplitude; "\
            "B_defect: resolution-decay of a variance/covariance "\
            "error-correlation effect) -- this restraint is a "\
            "physically motivated prior, not an identity. <= 0 (or "\
            "None) disables it."
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

def a_k_smoothness_penalty_and_gradient(theta, weight):
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
  (a_k_smoothness_weight); <= 0 or None returns a no-op.

  Returns (penalty, d(penalty)/d(theta)), penalty a plain float, the
  gradient a numpy array the same shape as theta (nonzero only in the
  a_k slots).
  """
  theta = np.asarray(theta, dtype=float)
  grad = np.zeros_like(theta)
  if(weight is None or weight <= 0):
    return 0.0, grad
  a, _, _ = dmodel.unpack_theta(theta)
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
  grad[:k] = g
  return penalty, grad

def b_defect_restraint_penalty_and_gradient(theta, b_sol_anchor, sigma):
  """ One-directional log-space restraint pulling B_defect toward
  b_sol_anchor (design doc sec. 6.4): penalty =
  0.5*(ln(B_defect) - ln(b_sol_anchor))^2 / sigma^2. Only B_defect's
  gradient component is nonzero -- b_sol_anchor is a fixed external
  value (the bulk-solvent step's own already-converged estimate, NOT a
  parameter being co-refined here), so there is no reverse pull on it.

  theta: natural-space parameter vector. b_sol_anchor: fixed positive
  float (the bulk-solvent fit's current B_sol point estimate). sigma:
  restraint width in ln(B) units (llgi_e_dmodel_params.
  b_sol_restraint_sigma); <= 0 or b_sol_anchor <= 0 returns a no-op.

  Returns (penalty, d(penalty)/d(theta)), penalty a plain float, the
  gradient a numpy array the same shape as theta (nonzero only in the
  B_defect slot, i.e. index -1).
  """
  theta = np.asarray(theta, dtype=float)
  grad = np.zeros_like(theta)
  if(sigma is None or sigma <= 0 or b_sol_anchor is None
     or b_sol_anchor <= 0):
    return 0.0, grad
  _, _, b_defect = dmodel.unpack_theta(theta)
  if(b_defect <= 0.0):
    return 0.0, grad
  log_diff = np.log(b_defect) - np.log(b_sol_anchor)
  penalty = 0.5 * (log_diff / sigma) ** 2
  grad[-1] = log_diff / (sigma * sigma * b_defect)
  return float(penalty), grad

class d_model_target_evaluator(object):
  """ scitbx.lbfgs target evaluator fitting D_model(s; theta) against
  the E-scale LLGI target, summed over the R-free/test set only (same
  restriction as mmtbx.refinement.llgi_e_bulk_solvent.
  e_sigmaa_target_evaluator), with Emodel (hence the bulk-solvent
  model) held fixed. Optimizes in unconstrained q-space (theta =
  exp(q), llgi_e_dmodel_reparam.py); .x is q, NOT theta directly.

  Uses scitbx.lbfgs's diag_mode="once" hook (design doc sec. 6.4): the
  diagonal of the q-space Hessian at the STARTING point seeds L-BFGS's
  initial inverse-Hessian, instead of the default isotropic guess --
  helps with theta's unevenly scaled natural parameters (a_k ~O(1) vs.
  B_defect ~O(10-300) in s^2 units; B_1..B_K are no longer part of
  theta at all -- a fixed ladder, b_k_grid, set at construction time --
  see llgi_e_dmodel.py's own module docstring for why). Entries that are
  not usefully positive at the start are replaced by a neutral 1 (see
  _diagonal_at).
  """

  # Diagonal curvatures at or below this are not usable for scaling the
  # first L-BFGS step (see _diagonal_at); scitbx.lbfgs's diag_mode hook
  # needs them strictly positive (tst_curvatures.py's
  # lbfgs_with_curvatures_mix_in.__call__).
  _curvature_floor = 1.e-6

  def __init__(self,
        e_eff, r_free_flags, e_model, dobs, centric_flags, d_star_sq,
        n_gaussian_terms=2, theta_start=None, max_iterations=200,
        b_sol_anchor=None,
        b_sol_restraint_sigma=0.6931471805599453,
        a_k_smoothness_weight=0.1, b_k_grid=None,
        include_constant_term=False):
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
    # LL/gradient/Hessian are normalized to a PER-REFLECTION MEAN below
    # (_natural_ll_grad_hess), matching llgi_e.h/ext.llgi_e_sigmaa_
    # target_and_gradients' own convention (target() is already divided
    # by n_selected) -- otherwise self.final_target/.target would be an
    # un-normalized sum over the test set, on a wildly different scale
    # from the spline path's reported target and not comparable to it
    # in logs/diagnostics (observed directly while testing this
    # integration: an un-normalized sum read as ~86000 next to the
    # spline path's ~-0.13 on the same data, even though the underlying
    # fits were equally well-behaved -- a presentation/consistency gap,
    # not a sign either fit was actually wrong).
    self.n_test = int(np.sum(test_sel))
    self.n_gaussian_terms = n_gaussian_terms
    self.b_sol_anchor = b_sol_anchor
    self.b_sol_restraint_sigma = b_sol_restraint_sigma
    self.a_k_smoothness_weight = a_k_smoothness_weight
    if(b_k_grid is None):
      b_k_grid = default_b_k_grid(n_gaussian_terms, self.s2)
      if(include_constant_term):
        b_k_grid = np.concatenate([[0.0], b_k_grid])
    self.b_k_grid = np.asarray(b_k_grid, dtype=float)
    n_terms = n_gaussian_terms + int(include_constant_term)
    assert self.b_k_grid.shape == (n_terms,)
    self.final_target = None

    if(theta_start is None):
      theta_start = self._default_theta_start(n_terms)
    theta_start = np.asarray(theta_start, dtype=float)
    assert theta_start.size == n_terms + 2
    self.x = flex.double(reparam.q_from_theta(theta_start))

    diag = self._diagonal_at(np.array(self.x))
    self.diag_mode = "once"
    self._diag_for_lbfgs = flex.double(1.0 / diag)

    term_parameters = scitbx.lbfgs.termination_parameters(
      max_iterations=max_iterations)
    exception_handling_parameters = scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound=True,
      ignore_line_search_failed_step_at_upper_bound=True)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=term_parameters,
      exception_handling_params=exception_handling_parameters)

  @staticmethod
  def _default_theta_start(k):
    """ Neutral starting theta: each a_k = 0.5/K (D_model itself is
    UNCONDITIONALLY bounded to (0,1) regardless of a_k -- see
    llgi_e_dmodel.py's own docstring -- so this is simply a modest,
    not-yet-committal starting scale, not a constraint-satisfying
    choice; B_k itself is no longer part of theta -- it is the fixed
    b_k_grid ladder, set separately, design doc sec. 6.4's well-
    posedness addendum), small initial b/B_defect (a modest, not-yet-
    committal defect term).
    """
    theta = np.empty(k + 2, dtype=float)
    if(k > 0):
      theta[:k] = 0.5 / k
    theta[-2] = 0.05
    theta[-1] = 40.0
    return theta

  def _natural_ll_grad_hess(self, q):
    theta = reparam.theta_from_q(q)
    ll, grad_p, hess_p = target.total_ll_gradient_hessian(
      theta, self.s2, self.e_eff, self.e_c, self.dobs,
      self.centric_flags, self.b_k_grid)
    # Normalize to a per-reflection mean (see __init__'s .n_test
    # comment) -- total_ll_gradient_hessian itself returns the raw,
    # un-normalized SUM (the natural quantity for Section 6's Fisher-
    # information work, where per-reflection curvature terms need to
    # add rather than be diluted), so the normalization is applied
    # here, at the evaluator/reporting boundary, not in that module.
    return theta, ll / self.n_test, grad_p / self.n_test, hess_p / self.n_test

  def _diagonal_at(self, q):
    theta, ll, grad_p, hess_p = self._natural_ll_grad_hess(q)
    diag_q = reparam.reparametrize_hessian_diagonal(
      theta, grad_p, np.diag(hess_p))
    # LL is being MAXIMIZED in natural/q space but scitbx.lbfgs
    # minimizes -- the diagonal fed to it must be the curvature of the
    # quantity actually being minimized, i.e. -LL, hence the sign flip
    # (see compute_functional_and_gradients_diag below for the same
    # flip applied to f/g).
    diag_q = -diag_q
    # Where the starting point is not locally convex in a parameter, use
    # a neutral unit curvature rather than the floor: 1/floor (1e6) as the
    # initial inverse-Hessian entry makes the first step so large that
    # the line search fails outright and the fit returns its starting
    # point (seen on 2G38 with k_mask=0 and no B_sol anchor).
    return np.where(diag_q > self._curvature_floor, diag_q, 1.0)

  def compute_functional_and_gradients(self):
    f, g, _ = self.compute_functional_gradients_diag()
    return f, g

  def compute_functional_gradients_diag(self):
    q = np.array(self.x)
    theta, ll, grad_p, hess_p = self._natural_ll_grad_hess(q)
    grad_q = reparam.reparametrize_gradient(theta, grad_p)

    b_penalty, b_grad_p = b_defect_restraint_penalty_and_gradient(
      theta, self.b_sol_anchor, self.b_sol_restraint_sigma)
    # Penalty gradient is in NATURAL (theta) space -- reparametrize via
    # the same dF/dq_i = dF/dp_i * p_i rule (llgi_e_dmodel_reparam.
    # reparametrize_gradient) before adding to grad_q, exactly as
    # grad_p itself was above.
    b_grad_q = reparam.reparametrize_gradient(theta, b_grad_p)

    smooth_penalty, smooth_grad_p = a_k_smoothness_penalty_and_gradient(
      theta, self.a_k_smoothness_weight)
    smooth_grad_q = reparam.reparametrize_gradient(theta, smooth_grad_p)

    # Minimize-me convention (matching every other evaluator in
    # llgi_e_bulk_solvent.py, and llgi_e.h's own target_one_h sign
    # flip): f = -(LL - penalties), penalties themselves already
    # minimize-me (added, not subtracted) by construction (see their
    # own docstrings). No negative-variance barrier or sum_k(a_k)<=1
    # constraint needed here (unlike an earlier version of this
    # evaluator): D_model(s; theta) is now UNCONDITIONALLY bounded to
    # (0, 1) by construction (tanh(smooth_relu(.))-wrapped, see
    # llgi_e_dmodel.py's own docstring for why a soft per-reflection
    # barrier could never actually fix the negative-variance issue --
    # the true log-likelihood diverges to -infinity approaching D=+-1,
    # so no finite barrier weight can outweigh the reward of crossing
    # into the negative-variance guard's region, which returns exactly
    # 0 -- and why sum_k(a_k)<=1 is now automatically guaranteed
    # regardless of a_k's own value, making that separate constraint
    # redundant too). The a_k smoothness penalty IS needed, unlike
    # those two: D_model=0 is a genuine, data-independent stationary
    # point of the raw likelihood (l(D=0)=0, l'(D=0)=0 always --
    # llgi_e_likelihood.py). With the OLDER free-B_k parametrization
    # that meant a_k->0/B_k->infinity was a real, reachable degenerate
    # direction (confirmed via tst_llgi_e_dmodel_fit.py); B_k is now a
    # fixed ladder (llgi_e_dmodel.py's own module docstring, design doc
    # sec. 6.4's well-posedness addendum), which already closes off
    # that specific runaway, but the fit is still free to send SOME
    # a_k individually toward 0 in a jagged, physically-implausible
    # pattern rather than genuinely fitting the data -- this penalty
    # discourages that in favour of a smooth coefficient profile across
    # the ladder, matching the coordinate-error-varies-smoothly-with-
    # B-factor physical picture directly.
    f = -ll + b_penalty + smooth_penalty
    g = -grad_q + b_grad_q + smooth_grad_q
    self.final_target = f

    if(getattr(self, "diag_mode", None) is not None):
      return f, flex.double(g), self._diag_for_lbfgs
    return f, flex.double(g)

  def theta(self):
    return reparam.theta_from_q(np.array(self.x))

def estimate_d_model_sigmaa(e_eff, r_free_flags, e_model, dobs,
      centric_flags, d_star_sq, n_gaussian_terms=2, max_iterations=200,
      b_sol_anchor=None,
      b_sol_restraint_sigma=0.6931471805599453,
      a_k_smoothness_weight=0.1, theta_start=None, b_k_grid=None,
      include_constant_term=False):
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
    b_sol_restraint_sigma=b_sol_restraint_sigma,
    a_k_smoothness_weight=a_k_smoothness_weight,
    b_k_grid=b_k_grid,
    include_constant_term=include_constant_term)
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
