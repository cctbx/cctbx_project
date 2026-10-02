from __future__ import absolute_import, division, print_function
import numpy as np

""" Unconstrained (log/exp) reparametrization of D_model(s; theta)'s
natural parameters theta = [a_1, B_1, ..., a_K, B_K, b, B_defect] (all
required positive -- sigmaA_model_handoff.md sec. 2), for use by an
unconstrained optimizer (L-BFGS). Every natural parameter here uses the
SAME p = exp(q) transform (the document's sec. 2 recommendation is
log-space for B_k/B_defect and an unspecified positive transform for
a_k/b "e.g. a_k = exp(alpha_k)" -- exp is used uniformly for all of
theta here, since every entry is positive-constrained the same way; the
soft sum_k a_k <= 1 constraint from sec. 2 is a separate penalty term,
not part of this reparametrization, and is not implemented in this
module).

Only the GRADIENT reparametrization is as simple as the document
suggests ("chain rule adjustments... are straightforward"): for p =
exp(q), dF/dq = dF/dp * p. The HESSIAN reparametrization is NOT simply
d2F/dp2 * p^2 -- there is an extra term:

    d2F/dq_i dq_j = d2F/dp_i dp_j * p_i * p_j  (+ dF/dp_i * p_i, i==j
                                                  only, diagonal term)

Full derivation for a multivariate p_i = exp(q_i) (each q_i
independent):

    dF/dq_i = dF/dp_i * p_i
    d2F/dq_i dq_j = d/dq_j [dF/dp_i * p_i]
                  = p_i * d/dq_j[dF/dp_i] + dF/dp_i * d(p_i)/dq_j
                  = p_i * p_j * d2F/dp_i dp_j + dF/dp_i * p_i * [i==j]

(the last term only survives on the diagonal, since d(p_i)/dq_j = 0 for
i != j -- p_i depends only on its own q_i). This extra diagonal term
only vanishes when the ORIGINAL-parameter gradient dF/dp_i is exactly
zero -- NOT guaranteed at an unconstrained (q-space) optimum for
parameters the design doc sec. 6.4 / sigmaA_model_handoff.md sec. 5
already expects to be poorly determined (near-flat, not exactly flat,
in the high-resolution B_k tail). Dropping this term was verified
(doc sec. 6.4) to be wrong by 3x or more away from convergence, and
still ~10% wrong fairly close to it -- do not drop it. See
llgi_e_dmodel_target.py for the natural-coordinate (p-space) gradient/
Hessian this module transforms.
"""

def q_from_theta(theta):
  """ q = ln(theta), elementwise. theta must be strictly positive
  (undefined/-inf otherwise -- callers are responsible for starting
  from a valid positive theta, e.g. phil defaults or a previous fit's
  converged values).
  """
  return np.log(np.asarray(theta, dtype=float))

def theta_from_q(q):
  """ theta = exp(q), elementwise. Always strictly positive by
  construction, for any finite real q -- the whole point of this
  reparametrization.

  q is clipped to +-_Q_CLIP before exponentiating -- an unconstrained
  L-BFGS line search can genuinely visit very large |q| while probing
  step lengths (observed in practice), and exp() of an unclipped large
  q silently overflows to inf, which then poisons every downstream
  D_model/likelihood computation with NaN.

  _Q_CLIP=30 (theta up to ~1e13) rather than a value nearer float64's
  true overflow boundary (exp(700)~1e304, exp(710) already overflows):
  a first version of this clip used 700, chosen only to avoid LITERAL
  exp() overflow, but that leaves theta itself free to reach ~1e304 --
  and theta this large overflows plenty of DOWNSTREAM consumers that
  multiply by theta directly (e.g. llgi_e_dmodel_reparam.
  reparametrize_gradient's grad_p*theta: even a merely tiny, not
  literally-zero grad_p, ~1e-300, still overflows against a ~1e304
  theta), confirmed on real 2G38 data as a second, independent source
  of overflow warnings after the D_model-side tanh-wrapping fix (see
  doc/llgi_target_design.md sec. 6.4) closed the first. 1e13 remains
  many orders of magnitude beyond any physically plausible a_k/B_k/b/
  B_defect value (design doc sec. 2's realistic ranges are O(1) to
  O(300)), so this still never constrains a sensible fit -- it only
  narrows the window in which a byproduct of clipping q can itself
  become numerically dangerous to whatever multiplies it next.
  """
  _Q_CLIP = 30.0
  q = np.clip(np.asarray(q, dtype=float), -_Q_CLIP, _Q_CLIP)
  return np.exp(q)

def reparametrize_gradient(theta, grad_p):
  """ dF/dq_i = dF/dp_i * p_i, elementwise. theta is the natural-space
  point (p, i.e. theta_from_q(q)) grad_p was evaluated at; grad_p is
  dF/dtheta (natural-space gradient, e.g. from llgi_e_dmodel_target.
  total_ll_gradient_hessian). Returns the q-space gradient, same shape.
  """
  theta = np.asarray(theta, dtype=float)
  grad_p = np.asarray(grad_p, dtype=float)
  return grad_p * theta

def reparametrize_hessian(theta, grad_p, hess_p):
  """ Full q-space Hessian, INCLUDING the extra diagonal term this
  module's own docstring derives (do not use hess_p * outer(theta,
  theta) alone -- see above). theta, grad_p as in
  reparametrize_gradient; hess_p is d2F/dtheta_i dtheta_j (natural-
  space Hessian, e.g. llgi_e_dmodel_target.total_ll_gradient_hessian's
  own hess return value). Returns the q-space Hessian, same shape as
  hess_p.
  """
  theta = np.asarray(theta, dtype=float)
  grad_p = np.asarray(grad_p, dtype=float)
  hess_p = np.asarray(hess_p, dtype=float)
  outer = np.outer(theta, theta)
  hess_q = hess_p * outer
  # Extra diagonal-only term: dF/dp_i * p_i, added to hess_q[i, i].
  diag_extra = grad_p * theta
  hess_q = hess_q + np.diag(diag_extra)
  return hess_q

def reparametrize_hessian_diagonal(theta, grad_p, hess_p_diagonal):
  """ Just the DIAGONAL of the q-space Hessian -- cheaper than the full
  reparametrize_hessian when only the diagonal is wanted (e.g. for
  scitbx.lbfgs's diag_mode preconditioning hook, design doc sec. 6.4).
  hess_p_diagonal is the natural-space Hessian's own diagonal
  (hess_p[i, i] for each i), not the full matrix. Returns a 1D array,
  same length as theta.
  """
  theta = np.asarray(theta, dtype=float)
  grad_p = np.asarray(grad_p, dtype=float)
  hess_p_diagonal = np.asarray(hess_p_diagonal, dtype=float)
  return hess_p_diagonal * theta * theta + grad_p * theta
