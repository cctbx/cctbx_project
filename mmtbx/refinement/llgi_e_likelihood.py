from __future__ import absolute_import, division, print_function
import numpy as np
from scipy.special import i0

""" Corrected per-reflection E-scale LLGI acentric/centric log-
likelihood and its first derivative w.r.t. D, for the D_model(s; theta)
sigmaA fit (see doc/llgi_target_design.md sec. 6.4).

These are the SYMMETRIC forms already used by cctbx_project/cctbx/xray/
targets/llgi_e.h's target_one_h/d_target_one_h_over_sigmaa (that C++ is
the reference implementation this module's l()/l_prime() must agree
with) -- NOT the forms literally written in sigmaA_model_handoff.md
sections 3.3/3.4/4.1/4.2, which were derived from the large-argument
ASYMPTOTIC approximation of ln(I_0(x)) (acentric) / log(cosh(x/2))
(centric) rather than the exact forms, and so do not match llgi_e.h
once the true, numerically stable Bessel/log-cosh functions are used
(see doc sec. 6.4 for the full derivation and finite-difference
verification of the correction).

Cross-checked directly against ext.llgi_e_sigmaa_target_and_gradients
(the real C++ functor) for a grid of (dobs, sigmaa, e_eff, e_model):
centric matches to machine precision; acentric TARGET matches to
~1e-7 rather than machine precision, traced to scitbx::math::bessel::
ln_of_i0's use of the classical Abramowitz & Stegun 9.8.1 rational-
polynomial approximation to I_0(x) (accurate to ~1.6e-7 by that
approximation's own known error bound) rather than an exact
evaluation -- acentric_l below uses scipy.special.i0 instead (via
np.log(i0(X))), which is essentially exact, for that piece; the
residual gap is a difference in NUMERICAL APPROXIMATION QUALITY of
ln(I_0(x)), not a disagreement about which formula is correct.
Acentric GRADIENT matches to MACHINE precision: _r(x) below (the Rice
ratio R(x)=I_1(x)/I_0(x) that acentric_l_prime depends on) is a direct port of scitbx::math::bessel::i1_over_i0's own
rational-polynomial approximation -- the exact same one llgi_e.h's own
gradient function uses -- not scipy.special.i1(x)/i0(x) (an earlier
version of this module used the latter; found to overflow to inf/nan
for x >~ 710 on real 2G38 data, since scipy's i0/i1 are each computed
as literal, individually huge numbers before dividing, unlike
i1_over_i0's two-separate-polynomials-then-divide approach, which never
forms an intermediate value that can overflow at any argument -- see
doc/llgi_target_design.md sec. 6.4).

Pure numpy/scipy functions only, no LLGI-target machinery, no theta/
D_model chain rule -- that combination is llgi_e_dmodel_target. This module operates entirely in terms of D (the
combined correlation coefficient D_c(h) = D_obs(h)*D_model(s_h; theta)
from sec. 3.2 of the handoff document), E_eff, E_c.

Negative-variance guard: llgi_e.h's target_one_h returns 0 (no
contribution) whenever v = 1-D^2 <= 0, rather than letting log(v) or
1/v blow up. D_model is bounded below SIGMAA_MAX, so this guard is a
safety net rather than a path the fit is expected to take. Every
function below applies the SAME s<=0 -> 0 guard as
llgi_e.h, so a caller summing over many reflections (llgi_e_dmodel_
target.py) gets a well-defined, finite (if temporarily useless)
contribution from an out-of-range reflection instead of a NaN that
poisons the whole sum -- matching llgi_e.h's own documented behaviour
in the one regime this module didn't originally handle.
"""

def _mask_invalid(s, result):
  """ Zero out `result` wherever s<=0 (matching llgi_e.h's v<=0 guard)
  OR result is itself non-finite for any other reason (defensive: e.g.
  extreme X overflowing i0(X) before s itself goes negative). Shared by
  every l/l_prime function below.
  """
  s = np.asarray(s, dtype=float)
  bad = (s <= 0.0) | ~np.isfinite(result)
  return np.where(bad, 0.0, result)

def _r(x):
  """ R(x) = I_1(x)/I_0(x) (the Rice ratio / Sim weighting function),
  computed as scitbx::math::bessel::i1_over_i0 does (cctbx_project/
  scitbx/math/bessel.h) -- via two SEPARATE rational-polynomial
  approximations to I_0(x) and I_1(x) (Abramowitz & Stegun 9.8.1/9.8.3:
  a small-x power series for |x|<3.75, a large-x asymptotic series
  otherwise), dividing only at the very end. Neither polynomial ever
  evaluates I_0(x)/I_1(x) as a literal huge number the way naive
  scipy.special.i0(x)/i1(x) does for large x (which genuinely overflows
  to inf past x~710, giving inf/inf=nan) -- this is the DIRECT-ratio
  approach, not an asymptotic patch bolted onto an overflowing
  computation (an earlier version of this function tried exactly that
  patch; ported the real function instead once it was noticed that
  scitbx already provides one, and that llgi_e.h's own gradient
  function already uses it for precisely this reason).
  """
  x = np.asarray(x, dtype=float)
  abs_x = np.abs(x)
  small = abs_x < 3.75
  # Small-x branch (Horner's-rule evaluation of A&S 9.8.1/9.8.3):
  y_small = np.where(small, x / 3.75, 0.0)
  y_small = y_small * y_small
  p = [1.0, 3.5156229, 3.0899424, 1.2067292, 0.2659732, 0.360768e-1,
       0.45813e-2]
  pp = [0.5, 0.87890594, 0.51498869, 0.15084934, 0.2658733e-1,
        0.301532e-2, 0.32411e-3]
  be0_small = np.zeros_like(x)
  be1_small = np.zeros_like(x)
  pow_y = np.ones_like(x)
  for i in range(7):
    be0_small = be0_small + p[i] * pow_y
    be1_small = be1_small + x * pp[i] * pow_y
    pow_y = pow_y * y_small
  # Large-x branch (A&S 9.8.2/9.8.4, in inverse powers of x):
  abs_x_safe = np.where(small, 1.0, abs_x)  # dummy where unused, discarded
  y_large = 3.75 / abs_x_safe
  q = [0.39894228, 0.1328592e-1, 0.225319e-2, -0.157565e-2, 0.916281e-2,
       -0.2057706e-1, 0.2635537e-1, -0.1647633e-1, 0.392377e-2]
  qq = [0.39894228, -0.3988024e-1, -0.362018e-2, 0.163801e-2,
        -0.1031555e-1, 0.2282967e-1, -0.2895312e-1, 0.1787654e-1,
        -0.420059e-2]
  be0_large = np.zeros_like(x)
  be1_large = np.zeros_like(x)
  pow_y = np.ones_like(x)
  for i in range(9):
    be0_large = be0_large + q[i] * pow_y
    be1_large = be1_large + qq[i] * pow_y
    pow_y = pow_y * y_large
  be0 = np.where(small, be0_small, be0_large)
  be1 = np.where(small, be1_small, be1_large)
  result = be1 / be0
  return np.where((x < 0.0) & (result > 0.0), -result, result)

def acentric_l(D, e_eff, e_c):
  """ l(D) = -ln(s) - D^2*(e_eff^2+e_c^2)/s + ln(I_0(X)),
  s = 1-D^2, X = 2*D*e_eff*e_c/s. The symmetric form matching llgi_e.h
  target_one_h's acentric branch (up to that function's own sign flip
  and dobs<=0/eeff<=0/emodel<=0 early-return guards, not reproduced
  here -- this module assumes valid, positive, physically meaningful
  inputs; callers are responsible for any such filtering, exactly as
  llgi_e.h itself documents).
  """
  D = np.asarray(D, dtype=float)
  s = 1.0 - D*D
  s_safe = np.where(s > 0.0, s, 1.0)  # dummy value where invalid, masked below
  quad = e_eff*e_eff + e_c*e_c
  X = 2.0*D*e_eff*e_c/s_safe
  result = -np.log(s_safe) - D*D*quad/s_safe + np.log(i0(X))
  return _mask_invalid(s, result)

def acentric_l_prime(D, e_eff, e_c):
  """ d(l)/dD, corrected symmetric-form closed formula (design doc sec.
  6.4): l'(D) = 2D/s - 2D*quad/s^2 + R(X)*k*(1+D^2)/s^2, k=2*e_eff*e_c.
  """
  D = np.asarray(D, dtype=float)
  s = 1.0 - D*D
  s_safe = np.where(s > 0.0, s, 1.0)
  quad = e_eff*e_eff + e_c*e_c
  k = 2.0*e_eff*e_c
  X = k*D/s_safe
  result = (2.0*D/s_safe - 2.0*D*quad/(s_safe*s_safe)
            + _r(X)*k*(1.0+D*D)/(s_safe*s_safe))
  return _mask_invalid(s, result)

def _log_cosh(y):
  """ Numerically stable log(cosh(y)) for y >= 0, matching llgi_e.h's
  centric branch: y + log((1+exp(-2y))/2), avoiding overflow in
  exp(y)/2 for large y (cosh(y) computed directly would overflow long
  before log(cosh(y)) itself becomes representable).
  """
  y = np.asarray(y, dtype=float)
  return y + np.log((1.0 + np.exp(-2.0*y)) / 2.0)

def centric_l(D, e_eff, e_c):
  """ l_c(D) = -0.5*[ln(s) + D^2*quad/s] + log_cosh(X/2), matching
  llgi_e.h's centric branch (ll_core/2 + ln_cosh(x/2)) -- the exact
  form, not sigmaA_model_handoff.md sec. 3.4's large-x-asymptotic
  Gaussian approximation (design doc sec. 6.4).
  """
  D = np.asarray(D, dtype=float)
  s = 1.0 - D*D
  s_safe = np.where(s > 0.0, s, 1.0)
  quad = e_eff*e_eff + e_c*e_c
  X = 2.0*D*e_eff*e_c/s_safe
  ll_core = -(np.log(s_safe) + D*D*quad/s_safe)
  result = ll_core/2.0 + _log_cosh(X/2.0)
  return _mask_invalid(s, result)

def centric_l_prime(D, e_eff, e_c):
  """ d(l_c)/dD, corrected symmetric-form closed formula (design doc
  sec. 6.4): l_c'(D) = D/s - D*quad/s^2 + 0.5*tanh(X/2)*dX/dD,
  dX/dD = k*(1+D^2)/s^2, k=2*e_eff*e_c.
  """
  D = np.asarray(D, dtype=float)
  s = 1.0 - D*D
  s_safe = np.where(s > 0.0, s, 1.0)
  quad = e_eff*e_eff + e_c*e_c
  k = 2.0*e_eff*e_c
  X = k*D/s_safe
  dXdD = k*(1.0+D*D)/(s_safe*s_safe)
  result = D/s_safe - D*quad/(s_safe*s_safe) + 0.5*np.tanh(X/2.0)*dXdD
  return _mask_invalid(s, result)
