from __future__ import absolute_import, division, print_function
import numpy as np

""" Physically-motivated D_model(s; theta) parametrization for the
E-scale LLGI sigmaA curve, replacing the B-spline-over-sigmoid fit
(mmtbx.refinement.llgi_e_bulk_solvent.e_sigmaa_target_evaluator) with a
closed form built from K positive Gaussian-decay terms (a discretized
coordinate-error B-factor distribution) minus one negative Gaussian term
(a bulk-solvent-mask-boundary-placement-error-induced negative
covariance between the ordered-atom and solvent-mask structure-factor
errors). See doc/llgi_target_design.md sec. 6.4 for the full design
discussion and its relationship to sigmaA_model_handoff.md sec. 2/4.3.

This module implements ONLY D_model(s; theta) and its first/second
partial derivatives w.r.t. theta (Section 7, implementation step 1) --
pure functions, no LLGI likelihood, no optimizer, no reparametrization.
Downstream steps (per-reflection likelihood/gradient/Hessian, chain-rule
combination, L-BFGS fit) build on top of this in later modules.

Parameter layout: theta is a flat sequence
    [a_1, a_2, ..., a_K, b, B_defect]
(K + 2 values total). B_1..B_K -- the K coordinate-error terms' decay
constants -- are NOT part of theta: they are a FIXED, externally
supplied ladder (b_k_grid, a plain array, passed as its own argument to
every function below), log-spaced across a physically plausible range
and never touched by the optimizer. Only the amplitudes a_k (linear in
D_raw) plus the defect term's own (b, B_defect) pair are fitted.

This is a deliberate reformulation, not the original design (design doc
sec. 6.4's "well-posedness" addendum) -- B_k started out as a fitted
parameter alongside a_k, matching sigmaA_model_handoff.md sec. 2's
original K-Gaussian-mixture proposal, but that is a textbook ILL-POSED
nonlinear fit: D_raw is a sum of radial basis functions all centered at
s2=0 with only their WIDTHS free, so (a) B_k enters D_raw nonlinearly
(dD_raw/dB_k = -a_k*s2*exp(-B_k*s2), a product of a_k and s2, unlike
a_k's own purely-linear dD_raw/da_k = exp(-B_k*s2)), and (b) two terms
with similar B_k are nearly linearly DEPENDENT (their exp(-B_k*s2)
curves are nearly identical shapes over any bounded s2 range), so nudge
one term's a_k up and the other's down without changing the sum at all
-- a genuinely flat direction, not merely a slowly-converging one. Both
together open a real, confirmed failure mode: since D_model(s;theta)->0
uniformly is a global stationary point of the raw likelihood for ANY
theta (l(D=0)=0 and l'(D=0)=0 identically -- llgi_e_likelihood.py), an
unrestrained free-B_k fit can and does run away along the (a_k->0,
B_k->infinity) direction -- an amplitude-vanishing, width-vanishing
Gaussian spike the fit becomes blind to, confirmed reproducibly on
tst_llgi_e_dmodel_fit.py's own synthetic data (B_k converging to
~1e13, hitting llgi_e_dmodel_reparam.theta_from_q's own numerical
_Q_CLIP ceiling, not a meaningful physical value) once D_model's high-
resolution asymptote was corrected to genuinely reach 0 (the fix
above); with the OLDER, buggy 0.5-asymptote wrapping, this direction
happened to be unreachable, accidentally masking the underlying
ill-posedness rather than fixing it.

The standard remedy for exactly this class of problem (fixed-center
RBF regression) is used here: STOP fitting the nonlinear B_k at all.
Put them on a fixed ladder spanning the physically plausible coordinate
-error range (chosen once, by the caller -- see
llgi_e_dmodel_fit.default_b_k_grid), and fit only the a_k amplitudes,
which enter D_raw LINEARLY -- this turns an ill-conditioned nonlinear
problem into a well-posed one (still needing a_k>=0 and, per the
physical picture of coordinate error varying smoothly across a
continuum of atomic B-factors, a mild smoothness prior across
neighbouring rungs of the ladder -- llgi_e_dmodel_fit.
a_k_smoothness_penalty_and_gradient), without reintroducing the
nonlinear-and-nearly-collinear B_k degrees of freedom that caused the
runaway in the first place. log/exp reparametrization for a_k/b/B_defect
is a separate concern, applied by the caller, not here (see doc sec.
6.4's note on why the Hessian reparametrization is NOT simply
multiplying by p^2 -- computing everything here in the natural
coordinates keeps that correction as a clean, separate, one-time step
rather than baked into these formulas).

D_model(s; theta) = tanh(smooth_relu(D_raw(s; theta))), where D_raw is
the plain sum of Gaussian terms (K+2 fitted params, as above, plus the
fixed b_k_grid ladder) and smooth_relu is a smooth approximation to
max(x, 0) (see _smooth_relu below). D_raw ITSELF is unbounded --
nothing in Section 2's a_k>=0 constraint (or the soft sum_k(a_k)<=1
penalty) stops it from exceeding 1 at some resolution even when it
starts below 1 at
s=0, and during an unconstrained L-BFGS fit it routinely does,
transiently, while the optimizer probes step lengths. That matters far
more than it sounds: the TRUE log-likelihood l(D) (llgi_e_likelihood.
py) diverges to -infinity as D->+-1 from inside the valid region
(s=1-D^2>0), while the negative-variance guard (s<=0) returns exactly
0 -- a discontinuous JUMP from -infinity up to 0 exactly at the
boundary, not a plateau. No finite-weight soft penalty on D_raw can
outweigh that (verified directly: at D=0.999, a weight=10, margin=0.95
quadratic barrier costs only ~0.024, while crossing to s<=0 is worth
+76.9 in target value) -- confirmed as the actual cause of transient
overflow warnings on real 2G38 data (macrocycle 1 of doc sec. 6.4's
15-cycle test run). The FIX, not just a mitigation: wrap D_raw so
D_model itself is confined for EVERY theta, unconditionally -- no
separate constraint/penalty machinery is needed at all, and
llgi_e_dmodel_fit.py's negative-variance barrier penalty is
accordingly removed (see that module).

Range is (0, 1), not (-1, 1), AND D_model(s -> infinity) = 0, not 0.5:
two constraints, one wrapping function, found the hard way across two
separate corrections. D_obs (the LLGI measurement-error correlation
coefficient the caller multiplies D_model by, sigmaA_model_handoff.md
sec. 3.2) is itself always non-negative by its own statistical
definition, and there is no case of practical interest where the true
model-error correlation D_model should be negative either -- so
D=D_obs*D_model should never need to range over negative values,
matching how the existing B-spline sigmaA fit already behaves (sigmaa
bounded to (0.01,0.99) via its own sigmoid, never negative). Separately,
D_model must -> 0 (not some other constant) as s -> infinity: every
term in D_raw is a decaying exponential, so D_raw(s -> infinity) = 0
identically, for ANY theta -- this is the correct, expected high-
resolution behaviour (no correlation between observed and calculated
amplitudes once resolution exceeds what the model can explain), and
the wrapping function must send D_raw=0 to D_model=0.

A first attempt at the (0,1)-range fix used D_model = 0.5*(1+tanh(.))
-- this DOES satisfy the non-negativity requirement, but sends
D_raw=0 to D_model=0.5, NOT 0, so EVERY fitted curve plateaued at 0.5
at high resolution regardless of theta -- a real, structural defect,
not a fitting failure, only caught after publishing a comparison of
curves that all showed exactly this shape (the user noticed the
flat-at-0.5 tail looked wrong and asked directly whether 0.5 was some
kind of built-in degenerate value). Plain (unrescaled) tanh(D_raw) has
the right asymptote (tanh(0)=0) but can go negative if the defect term
b*exp(-B_defect*s^2) exceeds the positive coordinate-error terms at
some s -- nothing in the current a_k>=0/b>=0/B_defect>=0 constraints
prevents that. The fix used here keeps both properties:
smooth_relu(x) is 0 at x=0 and effectively max(x,0) elsewhere (a
smoothed hinge, not a hard clip, so d_model_gradient/d_model_hessian
stay well-defined everywhere), so tanh(smooth_relu(D_raw)) is 0 at
D_raw=0 (correct high-resolution asymptote), never negative even if
D_raw dips below 0 (a negative D_raw is treated as "no positive
signal", not as a scientifically meaningful negative correlation), and
saturates toward 1 as D_raw grows positively.

d_model_raw/d_model_raw_gradient/d_model_raw_hessian (below) are the
original, UNBOUNDED sum-of-Gaussians functions, kept as the inner
building block; d_model/d_model_gradient/d_model_hessian apply
tanh(smooth_relu(.)) and are the versions actually consumed everywhere
downstream.
"""

def unpack_theta(theta):
  """ Split a flat theta array into (a, b_val, b_defect_val), each a 1D
  numpy array/scalar: a = [a_1..a_K] (the K coordinate-error term
  AMPLITUDES -- the corresponding decay constants B_1..B_K are NOT part
  of theta; they are a fixed, externally-supplied ladder, b_k_grid,
  passed separately to every function below -- see this module's own
  docstring for why), b = scalar (defect-term amplitude), B_defect =
  scalar (defect-term decay constant, still fitted, unlike B_1..B_K).
  Raises ValueError if theta's length is not K+2 for some non-negative
  integer K.
  """
  theta = np.asarray(theta, dtype=float)
  if(theta.size < 2):
    raise ValueError(
      "unpack_theta: theta must have length >= 2 (K+2 for K>=0 "
      "coordinate-error amplitudes plus the mandatory b/B_defect "
      "pair); got length %d." % theta.size)
  a = theta[:-2]
  b = theta[-2]
  b_defect = theta[-1]
  return a, b, b_defect

def d_model_raw(s2, theta, b_k_grid):
  """ D_raw(s; theta) = sum_k a_k*exp(-B_k*s2) - b*exp(-B_defect*s2),
  evaluated at s2 = s^2 (an array of any shape, e.g. one value per
  reflection or per resolution bin). b_k_grid is the FIXED ladder of
  B_1..B_K decay constants (length K, matching a's length) -- see this
  module's own docstring for why B_k is fixed rather than fitted.
  Returns an array the same shape as s2. UNBOUNDED -- see this module's
  own docstring for why the actually-consumed D_model(s; theta) =
  tanh(smooth_relu(D_raw(s; theta))) instead (below).
  """
  s2 = np.asarray(s2, dtype=float)
  a, b, b_defect = unpack_theta(theta)
  b_k_grid = np.asarray(b_k_grid, dtype=float)
  assert b_k_grid.shape == a.shape
  total = np.zeros_like(s2, dtype=float)
  for a_k, B_k in zip(a, b_k_grid):
    total = total + a_k * np.exp(-B_k * s2)
  total = total - b * np.exp(-b_defect * s2)
  return total

def d_model_raw_gradient(s2, theta, b_k_grid):
  """ dD_raw/dtheta_i, evaluated at s2 (array), for every parameter in
  theta, in the SAME order as theta itself (b_k_grid as in
  d_model_raw -- fixed, not part of theta, so there is no dD_raw/dB_k
  component here at all). Returns an array of shape (theta.size,) +
  s2.shape (i.e. one gradient component per theta parameter, each
  itself an array over s2's shape).

  Closed forms (design doc sec. 6.4 / sigmaA_model_handoff.md sec. 4.3):
    dD_raw/da_k       = exp(-B_k*s2)   (B_k fixed, from b_k_grid)
    dD_raw/db          = -exp(-B_defect*s2)
    dD_raw/dB_defect   = b*s2*exp(-B_defect*s2)
  """
  s2 = np.asarray(s2, dtype=float)
  theta = np.asarray(theta, dtype=float)
  a, b, b_defect = unpack_theta(theta)
  b_k_grid = np.asarray(b_k_grid, dtype=float)
  assert b_k_grid.shape == a.shape
  k = a.size
  grad = np.zeros((theta.size,) + s2.shape, dtype=float)
  for j in range(k):
    B_k = b_k_grid[j]
    grad[j] = np.exp(-B_k * s2)
  e_defect = np.exp(-b_defect * s2)
  grad[-2] = -e_defect
  grad[-1] = b * s2 * e_defect
  return grad

def d_model_raw_hessian(s2, theta, b_k_grid):
  """ d2D_raw/dtheta_i dtheta_j, evaluated at s2 (array), for every
  (i, j) pair (b_k_grid as in d_model_raw). Returns an array of shape
  (theta.size, theta.size) + s2.shape. D_raw is now LINEAR in every
  remaining fitted parameter (each a_k enters as a_k*exp(-B_k*s2) with
  B_k FIXED -- unlike the earlier free-B_k form, there is no a_k*B_k
  product term left at all), so every a_k/a_j second derivative
  (including a_k with itself) is exactly zero; only the defect term's
  own (b, B_defect) 2x2 block is nonzero, same as before:
    d2D_raw/da_k da_j (any i, j, including i==j) = 0
    d2D_raw/db^2               = 0
    d2D_raw/db dB_defect       = s2*exp(-B_defect*s2)
    d2D_raw/dB_defect^2        = -b*s2^2*exp(-B_defect*s2)
  """
  s2 = np.asarray(s2, dtype=float)
  theta = np.asarray(theta, dtype=float)
  a, b, b_defect = unpack_theta(theta)
  hess = np.zeros((theta.size, theta.size) + s2.shape, dtype=float)
  e_defect = np.exp(-b_defect * s2)
  d_b_bdef = s2 * e_defect
  d_bdef_bdef = -b * s2 * s2 * e_defect
  hess[-2, -1] = d_b_bdef
  hess[-1, -2] = d_b_bdef
  hess[-1, -1] = d_bdef_bdef
  return hess


# tanh(x) saturates to EXACTLY 1.0 in float64 once |x| exceeds roughly
# 19.06 (well beyond any physically plausible D_raw, but genuinely
# reachable during an unconstrained L-BFGS line search's step-length
# probing -- observed as the residual failure mode of the tanh
# approach itself: a saturated D_model would put s=1-D_model^2 exactly
# at 0, hitting the same negative-variance guard this wrapping exists
# to avoid, just far further out and with the gradient toward it
# already vanishing (1-tanh(x)^2 -> 0), rather than the original
# adversarial-reward failure mode). Clipping the smooth_relu(D_raw)
# argument to tanh keeps D_model strictly inside (0, 1) with a
# comfortable float64 margin, closing this off explicitly rather than
# relying on it being merely astronomically unlikely.
_D_RAW_CLIP = 15.0

# smooth_relu(x) = 0.5*(x + sqrt(x^2+eps)) -- a smooth approximation to
# max(x, 0): smooth_relu(0) ~ 0.5*sqrt(eps) (near-exactly 0 for a small
# eps, unlike a hard max(x,0) which would introduce a non-differentiable
# kink at x=0 that an L-BFGS fit could land on exactly), smooth_relu(x)
# -> x for x >> 0, smooth_relu(x) -> 0 for x << 0. _SMOOTH_RELU_EPS
# chosen small enough that smooth_relu(0) ~ 7e-4 is negligible next to
# any physically meaningful D_model difference, while keeping the
# transition sharp enough that D_model still closely tracks a genuine
# max(D_raw,0)-then-tanh shape rather than a visibly softened one.
_SMOOTH_RELU_EPS = 1.e-6

def _smooth_relu(x):
  return 0.5 * (x + np.sqrt(x * x + _SMOOTH_RELU_EPS))

def _smooth_relu_prime(x):
  return 0.5 * (1.0 + x / np.sqrt(x * x + _SMOOTH_RELU_EPS))

def _smooth_relu_double_prime(x):
  return 0.5 * _SMOOTH_RELU_EPS / np.power(x * x + _SMOOTH_RELU_EPS, 1.5)

def d_model(s2, theta, b_k_grid):
  """ D_model(s; theta) = tanh(smooth_relu(D_raw(s; theta))) -- the
  actual, BOUNDED-to-(0,1) sigmaA curve consumed everywhere downstream,
  with D_model -> 0 as D_raw -> 0 (the correct high-resolution
  asymptote: every term in D_raw decays to 0 as s -> infinity, for ANY
  theta) and D_model never negative even where D_raw dips below 0 (see
  this module's own docstring for the full derivation, including why a
  simpler 0.5*(1+tanh(.)) rescale was tried first and found wrong).
  Returns an array the same shape as s2.
  """
  # draw is clipped BEFORE _smooth_relu, not just its output: for
  # |draw| far beyond _D_RAW_CLIP (genuinely reached for extreme theta
  # during L-BFGS line-search probing -- confirmed on real data),
  # _smooth_relu's own x*x term can overflow before tanh ever sees it.
  # Clipping draw first is safe (_smooth_relu is monotonic and
  # tanh(smooth_relu(x)) is already saturated for |x| > _D_RAW_CLIP
  # regardless of exactly how large x is beyond that).
  draw = np.clip(d_model_raw(s2, theta, b_k_grid), -_D_RAW_CLIP, _D_RAW_CLIP)
  return np.tanh(_smooth_relu(draw))

def d_model_gradient(s2, theta, b_k_grid):
  """ dD_model/dtheta_i = (1-tanh(u)^2) * smooth_relu'(D_raw) *
  dD_raw/dtheta_i, where u = smooth_relu(D_raw) -- the tanh chain rule
  applied through the extra smooth_relu layer. Returns an array of
  shape (theta.size,) + s2.shape, verified against finite differences
  (see tst_llgi_e_dmodel.py).

  Where u = smooth_relu(D_raw) exceeds _D_RAW_CLIP in magnitude (D_raw
  itself far outside any physically sane range -- smooth_relu(D_raw)
  and D_raw coincide for D_raw > 0, so this is the same "extreme theta
  during L-BFGS line-search probing" regime the original tanh(D_raw)
  wrapping had to guard against), the gradient is forced to exactly 0
  rather than computed from the formula above -- NOT an approximation:
  in the TRUE (unclipped) limit, sech2=1-tanh(u)^2 decays like
  exp(-2*|u|), which beats any of d_model_raw_gradient's own terms
  (each at most polynomial in theta for fixed s2), so the true product
  genuinely -> 0 as |u| grows (same reasoning verified for the earlier
  tanh(D_raw)-only wrapping; smooth_relu does not change this, since
  smooth_relu(x)=x for x>0 and is bounded near 0 for x<=0). But u gets
  CLIPPED for tanh's own sake, which freezes sech2 at a small-but-fixed
  value while d_model_raw_gradient (computed from the real, unclipped
  theta) can still grow without bound -- multiplying a frozen-small
  sech2 by an unbounded raw gradient does NOT reproduce the true (zero)
  limit, it diverges. Forcing the gradient to 0 exactly where clipping
  was triggered restores the correct limiting behaviour.
  """
  s2 = np.asarray(s2, dtype=float)
  draw_raw = d_model_raw(s2, theta, b_k_grid)
  # Test for "extreme" on draw_raw DIRECTLY (a cheap comparison, no
  # overflow risk) before ever calling _smooth_relu on it --
  # _smooth_relu's own x*x term can itself overflow for extreme
  # draw_raw (confirmed: 1e300-scale a_k on real 2G38 data), so
  # draw must be clipped BEFORE _smooth_relu sees it, not after.
  clipped = np.abs(draw_raw) > _D_RAW_CLIP
  draw = np.clip(draw_raw, -_D_RAW_CLIP, _D_RAW_CLIP)
  u = _smooth_relu(draw)
  t = np.tanh(u)
  sech2 = 1.0 - t * t
  srp = _smooth_relu_prime(draw)
  # d_model_raw_gradient can itself overflow to inf where clipped is
  # True (see this function's own docstring) -- the result is masked
  # to exactly 0 there regardless, so the overflow is provably benign,
  # not silently hidden; suppress only the warning it would otherwise
  # raise, not the (already-verified-correct) masking logic below.
  with np.errstate(over="ignore", invalid="ignore"):
    grad = ((sech2 * srp)[np.newaxis, ...]
            * d_model_raw_gradient(s2, theta, b_k_grid))
  return np.where(clipped[np.newaxis, ...], 0.0, grad)

def d_model_hessian(s2, theta, b_k_grid):
  """ d2D_model/dtheta_i dtheta_j, tanh chain rule applied through the
  extra smooth_relu layer (NOT simply (1-tanh(u)^2)*smooth_relu'(D_raw)
  times d_model_raw_hessian -- see this module's own docstring on why
  an analogous omission in the q=exp(theta) reparametrization elsewhere
  in this design was wrong by up to 3x-40%+ away from a stationary
  point; the same hazard applies here, TWICE over -- once for the tanh
  layer, once for the smooth_relu layer). Writing u = smooth_relu
  (D_raw):

    du/dtheta_i = smooth_relu'(D_raw) * dD_raw/dtheta_i
    d2u/dtheta_i dtheta_j = smooth_relu''(D_raw) * dD_raw/dtheta_i *
        dD_raw/dtheta_j + smooth_relu'(D_raw) * d2D_raw/dtheta_i dtheta_j
    d2D_model/dtheta_i dtheta_j = -2*tanh(u)*(1-tanh(u)^2)
        * du/dtheta_i * du/dtheta_j + (1-tanh(u)^2) * d2u/dtheta_i dtheta_j

  Returns an array of shape (theta.size, theta.size) + s2.shape,
  verified against finite differences (see tst_llgi_e_dmodel.py). Where
  |u| exceeds _D_RAW_CLIP, forced to exactly 0 -- see d_model_gradient's
  own docstring for why (the same true-limit-is-zero-but-clipping-
  breaks-the-cancellation reasoning applies here, for every term).
  """
  s2 = np.asarray(s2, dtype=float)
  draw_raw = d_model_raw(s2, theta, b_k_grid)
  # Same clip-before-_smooth_relu reasoning as d_model_gradient (see
  # its own docstring) -- test draw_raw directly (cheap, no overflow
  # risk) before ever calling _smooth_relu on it.
  clipped = np.abs(draw_raw) > _D_RAW_CLIP
  draw = np.clip(draw_raw, -_D_RAW_CLIP, _D_RAW_CLIP)
  u = _smooth_relu(draw)
  t = np.tanh(u)
  sech2 = 1.0 - t * t
  srp = _smooth_relu_prime(draw)
  srpp = _smooth_relu_double_prime(draw)
  # Same provably-benign-overflow situation as d_model_gradient (see
  # its own docstring) -- g/h can overflow where clipped is True, but
  # the result is masked to exactly 0 there regardless.
  with np.errstate(over="ignore", invalid="ignore"):
    g = d_model_raw_gradient(s2, theta, b_k_grid)  # shape (p,) + s2.shape
    h = d_model_raw_hessian(s2, theta, b_k_grid)   # shape (p, p) + s2.shape
    du = srp[np.newaxis, ...] * g           # shape (p,) + s2.shape
    d2u = (srpp[np.newaxis, np.newaxis, ...]
           * (g[:, np.newaxis, ...] * g[np.newaxis, :, ...])
           + srp[np.newaxis, np.newaxis, ...] * h)
    outer_du = du[:, np.newaxis, ...] * du[np.newaxis, :, ...]
    term1 = (-2.0 * t * sech2)[np.newaxis, np.newaxis, ...] * outer_du
    term2 = sech2[np.newaxis, np.newaxis, ...] * d2u
    hess = term1 + term2
  return np.where(clipped[np.newaxis, np.newaxis, ...], 0.0, hess)
