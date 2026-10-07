from __future__ import absolute_import, division, print_function
import numpy as np

""" D_model(s; theta): the parametrised sigmaA(resolution) curve fitted on
the E scale (mmtbx.refinement.llgi_e_dmodel_fit). With s2 = s^2 = d*^2/4,

  D_raw(s2)   = sum_k a_k*exp(-B_k*s2) - b*exp(-B_defect*s2)
  D_model(s2) = SIGMAA_MAX*tanh(smooth_relu(D_raw(s2)))

The positive terms stand for the distribution of coordinate errors, the
negative one for the error in the bulk-solvent boundary, whose structure
factor errors are anticorrelated with those of the atoms. theta =
[a_1..a_K, b, B_defect]. The B_k are a fixed ladder (b_k_grid, not part of
theta), so D_raw is linear in the a_k; fitting the B_k as well is ill-posed
(terms with similar B_k are nearly collinear). The ladder may include
B = 0, a resolution-independent term.

The wrapping keeps D_model in (0, SIGMAA_MAX) for any theta, so that
D = Dobs*D_model stays below 1 (the likelihood diverges as D -> 1), and
smooth_relu keeps it non-negative where the defect term dominates. Without
a B = 0 term, D_model -> 0 at high resolution.

Pure numpy functions: values and first derivatives with respect to theta.
"""

def unpack_theta(theta):
  """ Split theta into (a, b, B_defect), a = [a_1..a_K]. Raises ValueError
  if theta has fewer than 2 elements.
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
  """ D_raw at s2 (any shape), for the fixed ladder b_k_grid (length K).
  Unbounded.
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
  """ dD_raw/dtheta, shape (theta.size,) + s2.shape:
    dD_raw/da_k       = exp(-B_k*s2)
    dD_raw/db         = -exp(-B_defect*s2)
    dD_raw/dB_defect  = b*s2*exp(-B_defect*s2)
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


# D_raw is clipped to +-_D_RAW_CLIP before smooth_relu and tanh: tanh is
# already saturated there (tanh(15) = 1 - 2e-13), and extreme theta during
# line searches would otherwise overflow x*x in smooth_relu.
_D_RAW_CLIP = 15.0

# smooth_relu(x) = 0.5*(x + sqrt(x^2 + eps)), a smooth max(x, 0);
# smooth_relu(0) = 5e-4.
_SMOOTH_RELU_EPS = 1.e-6

def _smooth_relu(x):
  return 0.5 * (x + np.sqrt(x * x + _SMOOTH_RELU_EPS))

def _smooth_relu_prime(x):
  return 0.5 * (1.0 + x / np.sqrt(x * x + _SMOOTH_RELU_EPS))

# sigmaA is capped below 1: the Rice and exact likelihoods both degenerate
# as sigmaA -> 1 (variance 1 - sigmaA^2 -> 0), and the bias-reduced map
# coefficients divide by 1 - (1 - phi)*sigmaA^2. Applied as a smooth scale
# on the tanh, so D_model and its derivatives stay smooth. It also keeps
# sigmaA below the hybrid LLGI's exact/Rice switch at 0.999.
SIGMAA_MAX = 0.995

def d_model(s2, theta, b_k_grid):
  """ D_model at s2 (any shape), in (0, SIGMAA_MAX). """
  draw = np.clip(d_model_raw(s2, theta, b_k_grid), -_D_RAW_CLIP, _D_RAW_CLIP)
  return SIGMAA_MAX * np.tanh(_smooth_relu(draw))

def d_model_gradient(s2, theta, b_k_grid):
  """ dD_model/dtheta, shape (theta.size,) + s2.shape. Zero where D_raw
  is clipped (|D_raw| > 15): there the true derivative is below 1.2e-9
  times dD_raw/dtheta, and dD_raw/dtheta itself can overflow.
  """
  s2 = np.asarray(s2, dtype=float)
  draw_raw = d_model_raw(s2, theta, b_k_grid)
  clipped = np.abs(draw_raw) > _D_RAW_CLIP
  draw = np.clip(draw_raw, -_D_RAW_CLIP, _D_RAW_CLIP)
  u = _smooth_relu(draw)
  t = np.tanh(u)
  sech2 = 1.0 - t * t
  srp = _smooth_relu_prime(draw)
  # overflow where clipped is masked below
  with np.errstate(over="ignore", invalid="ignore"):
    grad = ((sech2 * srp)[np.newaxis, ...]
            * d_model_raw_gradient(s2, theta, b_k_grid))
  return SIGMAA_MAX * np.where(clipped[np.newaxis, ...], 0.0, grad)
