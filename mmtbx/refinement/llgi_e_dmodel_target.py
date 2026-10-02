from __future__ import absolute_import, division, print_function
import numpy as np
import mmtbx.refinement.llgi_e_dmodel as dmodel
import mmtbx.refinement.llgi_e_likelihood as lik

""" Chain-rule combination of D_model(s; theta) (llgi_e_dmodel.py) and
the per-reflection E-scale LLGI log-likelihood (llgi_e_likelihood.py)
into a full theta-space log-likelihood, gradient, and Hessian, summed
over a set of reflections (Section 7, implementation step 3 -- see
doc/llgi_target_design.md sec. 6.4).

D_c(h) = D_obs(h) * D_model(s_h; theta)  (sigmaA_model_handoff.md sec.
3.2) is the only place theta enters the likelihood. Per reflection:

    dLL_h/dtheta_i = l'(D_c) * D_obs(h) * g_i(s_h)
    d2LL_h/dtheta_i dtheta_j = l''(D_c) * D_obs(h)^2 * g_i(s_h)*g_j(s_h)
                                 + l'(D_c) * D_obs(h) * h_ij(s_h)

using l', l'' from llgi_e_likelihood's acentric or centric forms as
appropriate per reflection (sigmaA_model_handoff.md sec. 4.3). LL(theta)
itself is sum_h l(D_c(h)) -- the quantity to MAXIMIZE (un-negated
log-likelihood-gain convention, matching llgi_e_likelihood.py; the
sign flip to a minimize-me convention, if wanted for an optimizer, is
the caller's responsibility, exactly as llgi_e.h's own target functor
applies it at its own seam).

Pure numpy functions only; no optimizer, no reparametrization, no
restraints (curvature penalty, B_sol/B_defect restraint) -- those are a
separate, later piece (the eventual L-BFGS target evaluator class).
"""

def total_ll_gradient_hessian(theta, s2, e_eff, e_c, dobs, centric_flags,
      b_k_grid):
  """ Sum LL(theta), its gradient, and its Hessian over every reflection
  in the input arrays (all 1D, same length n_refl; s2 = s^2 per
  reflection, e_eff/e_c/dobs the fixed per-reflection LLGI inputs,
  centric_flags a boolean array selecting the centric formula per
  reflection). b_k_grid: the fixed ladder of B_1..B_K coordinate-error
  decay constants (llgi_e_dmodel.py's own module docstring) -- NOT part
  of theta, passed through unchanged to every dmodel.d_model* call.

  Returns (LL, grad, hess): LL a scalar, grad shape (theta.size,), hess
  shape (theta.size, theta.size).
  """
  s2 = np.asarray(s2, dtype=float)
  e_eff = np.asarray(e_eff, dtype=float)
  e_c = np.asarray(e_c, dtype=float)
  dobs = np.asarray(dobs, dtype=float)
  centric_flags = np.asarray(centric_flags, dtype=bool)
  theta = np.asarray(theta, dtype=float)
  n = s2.shape[0]
  assert e_eff.shape == (n,)
  assert e_c.shape == (n,)
  assert dobs.shape == (n,)
  assert centric_flags.shape == (n,)

  d_model_vals = dmodel.d_model(s2, theta, b_k_grid)  # shape (n,)
  g = dmodel.d_model_gradient(s2, theta, b_k_grid)     # shape (p, n)
  h = dmodel.d_model_hessian(s2, theta, b_k_grid)      # shape (p, p, n)
  D = dobs * d_model_vals                             # shape (n,)

  ll_per_refl = np.empty(n, dtype=float)
  lp_per_refl = np.empty(n, dtype=float)
  lpp_per_refl = np.empty(n, dtype=float)

  acentric_sel = ~centric_flags
  if(np.any(acentric_sel)):
    Da, ea, ca = D[acentric_sel], e_eff[acentric_sel], e_c[acentric_sel]
    ll_per_refl[acentric_sel] = lik.acentric_l(Da, ea, ca)
    lp_per_refl[acentric_sel] = lik.acentric_l_prime(Da, ea, ca)
    lpp_per_refl[acentric_sel] = lik.acentric_l_double_prime(Da, ea, ca)
  if(np.any(centric_flags)):
    Dc, ec_eff, ec_c = (
      D[centric_flags], e_eff[centric_flags], e_c[centric_flags])
    ll_per_refl[centric_flags] = lik.centric_l(Dc, ec_eff, ec_c)
    lp_per_refl[centric_flags] = lik.centric_l_prime(Dc, ec_eff, ec_c)
    lpp_per_refl[centric_flags] = lik.centric_l_double_prime(
      Dc, ec_eff, ec_c)

  LL = float(np.sum(ll_per_refl))

  # dLL_h/dtheta_i = l'(D)*dobs*g_i(s_h); sum over h.
  grad = np.einsum('n,n,pn->p', lp_per_refl, dobs, g)

  # d2LL_h/dtheta_i dtheta_j = l''(D)*dobs^2*g_i*g_j + l'(D)*dobs*h_ij;
  # sum over h.
  hess = (
    np.einsum('n,n,pn,qn->pq', lpp_per_refl, dobs*dobs, g, g)
    + np.einsum('n,n,pqn->pq', lp_per_refl, dobs, h))

  return LL, grad, hess
