from __future__ import absolute_import, division, print_function
""" The SigmaN model of the LLGI data normalisation: the expected intensity
of the data (without the epsilon factor) as

  ln SigmaN(h) = ln_scale + ln BEST(s^2) + ln b_bin(h)
                 - h.beta.h - x(h) h.P(beta1).h,

with BEST the X-ray average-intensity curve of phasertng (best_curve_at_ssqr),
b_bin a per-resolution-bin correction to its shape, beta the anisotropy
tensor (in phasertng's convention, exp(-(b11 h^2 + b22 k^2 + b33 l^2 + b12 hk
+ b13 hl + b23 kl)), reduced to its independent parameters by the space
group's symmetry) and beta1 an optional resolution-dependent anisotropy,
multiplied by x = d*^2/d*^2_max and with its isotropic part removed by
P (SigmaN::AnisoBetaRemoveIsoB), so that the isotropic fall-off stays with
the BEST curve and the bin corrections.

The model has the restraints of phasertng's RefineANISO: a sphericity
restraint on the anisotropic part of beta (at the resolution limit, where
the resolution-dependent part is largest) and of beta1, and a restraint on
ln b_bin to zero, weak at low resolution (sigma = 0.15 + 0.001 d^3).
Parameters are [ln_scale, ln_b_bin[0..n-1], beta_ind[...], beta1_ind[...]].

fit() refines this model jointly with a shared tNCS coherence against the
intensities alone, with the Wilson likelihood including the measurement
error (cctbx.xray.ext.wilson_intensity): the expected intensity is
epsilon SigmaN TEPS, TEPS = max(1 + Gamma(s) sum_g A_g, TEPS_FLOOR) with
A_g the consensus pair terms of the model's tNCS groups and Gamma(s) =
gamma exp(-2 pi^2 v^2 s^2) (llgi_tncs.coherence_curve). The TEPS returned
is floored at U_MIN instead, as it is the frame used downstream.
"""
import math
import numpy
import iotbx.phil
from libtbx import group_args
from libtbx.utils import null_out
from cctbx import adptbx, sgtbx
from cctbx.array_family import flex
from cctbx.xray import ext as xray_ext
from mmtbx.refinement.llgi_tncs import coherence_curve, U_MIN

# numerical guard only, because the measurement error already removes the
# influence of reflections with Sigma << sigma_I
TEPS_FLOOR = 1e-4

# phasertng io/data/best_curve_at_ssqr.h (X-ray BEST table), verbatim:
# intensity at ssqr = BEST_SSQR0 + i*BEST_STEP, to be multiplied by BEST_K
BEST_SSQR0 = 0.009
BEST_STEP = 0.004091973
BEST_K = 8.563714792513845e-06
BEST_TABLE = (
  119420.094, 101659.641, 81690.414, 63365.090, 57875.355,
  61608.578, 72501.570, 90078.656, 109741.945, 127082.102,
  137412.297, 139580.219, 137333.297, 133472.188, 129507.602,
  125845.414, 120904.797, 115925.875, 110775.266, 104672.258,
  97298.055, 89829.648, 83173.445, 78112.539, 73416.969,
  69432.289, 65371.203, 61651.734, 58562.309, 55816.238,
  53610.000, 52034.773, 51116.430, 50133.047, 49400.891,
  48117.367, 47040.086, 46641.746, 46537.227, 45864.223,
  46024.473, 45880.148, 46072.418, 45931.863, 46082.207,
  46202.676, 46375.520, 45861.191, 45652.387, 45039.852,
  44479.777, 43549.965, 43501.113, 43396.234, 42620.055,
  41100.250, 39922.781, 39130.934, 38021.281, 36782.801,
  35499.730, 34165.281, 32979.105, 31575.912, 30199.180,
  29096.426, 28114.213, 27363.445, 26257.604, 25265.410,
  24337.518, 23700.688, 22804.635, 22165.574, 21679.498,
  21230.615, 20727.221, 20303.465, 20050.027, 19911.986,
  19610.537, 19203.639, 18901.963, 18540.625, 18118.111,
  17952.115, 17938.740, 17940.045, 17767.729, 17619.865,
  17590.848, 17258.307, 17119.484, 16943.379, 16874.393,
  16605.168, 16426.338, 16412.588, 16348.048, 16272.690,
  15996.186, 15787.716, 15585.001, 15554.448, 15461.621,
  15466.696, 15320.923, 15018.120, 14736.160, 14501.164,
  14441.944, 14253.150, 14107.798, 13894.836, 13780.529,
  13684.477, 13461.170, 13299.096, 13077.154, 13108.808,
  13084.907, 13235.410, 13108.882, 13132.099, 13024.892,
  12993.081, 12830.820, 12948.324, 12934.137, 12931.499,
  12821.528, 12988.473, 13004.402, 13063.121, 13034.793,
  13068.688, 13004.556, 12888.828, 12963.618, 13111.302,
  13113.813, 13059.174, 13050.271, 13214.941, 13377.753,
  13383.613, 13239.169, 13218.517, 13199.132, 13154.157,
  13054.874, 13083.780, 13230.734, 13353.612, 13294.563,
  13324.932, 13357.476, 13570.960, 13736.209, 13983.591,
  14122.543, 14225.675, 14095.190, 14082.426, 14167.668,
  14348.723, 14491.771, 14480.615, 14510.075, 14619.429,
  14729.702, 14846.838, 14874.652, 14945.954, 15035.918,
  15075.190, 15067.756, 14902.410, 14830.787, 14857.209,
  14914.312, 14954.278, 15112.135, 15209.297, 15205.097,
  15024.275, 14929.255, 14948.889, 14968.661, 14960.080,
  14879.511, 14852.724, 14716.297, 14628.658, 14479.229,
  14316.959, 14116.013, 14189.159, 14226.442, 13770.954,
  13480.519, 13410.201, 13314.349, 13192.502, 13069.858,
  13146.034, 13034.913, 12845.740, 12565.674, 12494.399,
  12431.720, 12224.740, 12045.704, 11936.945, 12001.480,
  12094.969, 12074.801, 12003.505, 11495.983, 11344.455,
  11263.591, 11174.380, 11040.017, 10927.041, 10810.141,
  10751.963, 10636.446, 10556.304, 10401.776, 10351.462,
  10445.242, 10447.674, 10341.939, 10140.619, 10027.756,
  9964.771, 9847.444, 9818.661, 9778.909, 9726.898,
  9673.646, 9495.260, 9443.504, 9342.603, 9237.232,
  9147.548, 8949.440, 9064.902, 8936.296, 8842.546,
  8568.799, 8741.543, 8745.086, 8737.993, 8733.811,
  8556.928, 8571.162, 8452.438, 8358.330, 8383.448,
  8434.308, 8541.982, 8522.707, 8395.525, 8380.148,
  8368.124, 8375.411, 8058.358, 7890.074, 7953.468,
  8102.355, 8013.817, 7888.717, 7916.124, 7980.988,
  8043.755, 7990.078, 7948.800, 8006.895, 7844.749,
  7776.404, 7710.022, 7620.016, 7512.692, 7396.036,
  7360.436, 7381.391, 7391.329, 7450.565, 7304.477,
  7169.270, 6941.855, 6925.536, 6893.606, 7031.191,
  7002.648, 6981.484, 6899.487, 6920.931, 6939.107
)

# phasertng RefineANISO: LOG_WILSON_SIGMA_CONSTANT, LOG_WILSON_SIGMA_FACTOR
LOG_WILSON_SIGMA_CONSTANT = 0.15
LOG_WILSON_SIGMA_FACTOR = 0.001

_BEST_TABLE_NP = numpy.array(BEST_TABLE)
assert len(BEST_TABLE) == 300

def _as_numpy(a, dtype=float):
  """ numpy copy of a flex array (or sequence) """
  if hasattr(a, "as_numpy_array"): return a.as_numpy_array().astype(dtype)
  return numpy.asarray(a, dtype=dtype)

def best_curve_values(d_star_sq):
  """ BEST curve (times k) at d*^2, as phasertng's best_curve_at_ssqr::get:
  linear interpolation, truncating the table index towards zero, constant
  for an index below 0 and from 299 on. flex.double in, flex.double out """
  ssqr = _as_numpy(d_star_sq)
  s1 = numpy.trunc((ssqr - BEST_SSQR0) / BEST_STEP).astype(int)
  low = s1 < 0
  high = s1 >= 299
  i1 = numpy.clip(s1, 0, 298)
  ssqr1 = i1 * BEST_STEP + BEST_SSQR0
  ssqr2 = (i1 + 1) * BEST_STEP + BEST_SSQR0
  b1 = _BEST_TABLE_NP[i1]
  b2 = _BEST_TABLE_NP[i1 + 1]
  v = b1 + (ssqr - ssqr1) * (b2 - b1) / (ssqr2 - ssqr1)
  v[low] = BEST_TABLE[0]
  v[high] = BEST_TABLE[299]
  return flex.double(BEST_K * v)

def bin_restraint_sigmas(d_mid):
  """ sigma of the restraint of ln b_bin to zero, X-ray weighting of
  phasertng RefineANISO: 0.15 + 0.001 d_mid^3 (d_mid in A) """
  return flex.double([LOG_WILSON_SIGMA_CONSTANT +
    LOG_WILSON_SIGMA_FACTOR * d ** 3 for d in d_mid])

def resolution_bins(d_star_sq, width=500, n_min=6, n_max=50, min_per_bin=10):
  """ Resolution bins as phasertng Bin::setup (default CHOICE 2): equal
  width in s^2 = d*^2 between the limits of the data, n = max(n_min,
  min(N // width, n_max)) of them, reduced by one while any bin has fewer
  than min_per_bin reflections (down to a single bin). Returns group_args
  (index: flex.int per reflection, n_bins, d_mid: flex.double of the
  mid-resolution in A, 1/sqrt of the mean d*^2 of the bin edges) """
  ssqr = _as_numpy(d_star_sq)
  n_refl = ssqr.size
  assert n_refl > 0
  ppm = 1.0e-6
  s = numpy.sqrt(ssqr)
  s_max = s.max() + ppm
  s_min = s.min()
  assert s_max > s_min + ppm
  n = n_min if n_min == n_max else max(n_min, min(n_refl // width, n_max))
  shell = (s_max ** 2 - s_min ** 2) / n
  s_max = math.ceil(s_max / ppm) * ppm
  s_min = math.floor(s_min / ppm) * ppm
  fns_min = s_min ** 2
  while True:
    index = numpy.floor((ssqr - fns_min) / shell).astype(int)
    index[s < s_min] = 0
    index = numpy.clip(index, 0, n - 1)  # phasertng: n for s < s_max -> n-1
    counts = numpy.bincount(index, minlength=n)
    if n == 1 or counts.min() >= min_per_bin: break
    n -= 1
    shell = (s_max ** 2 - s_min ** 2) / n
  d_mid = 1 / numpy.sqrt(fns_min + (numpy.arange(n) + 0.5) * shell)
  return group_args(index=flex.int(index.astype(numpy.int32)), n_bins=n,
    d_mid=flex.double(d_mid))

# beta, in phasertng's order and convention: exp(-(b11 h^2 + b22 k^2 +
# b33 l^2 + b12 hk + b13 hl + b23 kl)); cctbx's tensors have the off-diagonal
# elements counted twice in the quadratic form
_TO_PHASERTNG = numpy.diag([1., 1., 1., 2., 2., 2.])

class sigman_model(object):
  """ SigmaN (without epsilon) as a linear function of the parameters
  x = [ln_scale, ln_b_bin[0..n-1], beta_ind, beta1_ind (if
  resolution_dependent)], plus ln of the BEST curve: ln_sigman(x) =
  d_ln_sigman_dx(x).x + ln BEST. Anisotropy (beta, beta1) is constrained by
  the symmetry of the space group """

  def __init__(self, miller_set, bins, resolution_dependent=False,
      resolution_dependent_restraint=10.):
    self.n_bins = bins.n_bins
    self.resolution_dependent = bool(resolution_dependent)
    self.resolution_dependent_restraint = float(resolution_dependent_restraint)
    self._bin_index = _as_numpy(bins.index, dtype=int)
    self._eps = _as_numpy(miller_set.epsilons().data())
    self._sigma_bin = _as_numpy(bin_restraint_sigmas(bins.d_mid))
    uc = miller_set.unit_cell()
    d_star_sq = miller_set.d_star_sq().data()
    ssqr = _as_numpy(d_star_sq)
    self._ln_best = numpy.log(_as_numpy(best_curve_values(d_star_sq)))
    ssqr_max = ssqr.max()
    self._x_res = ssqr / ssqr_max
    hkl = _as_numpy(miller_set.indices().as_vec3_double()).reshape(-1, 3)
    h, k, l = hkl[:, 0], hkl[:, 1], hkl[:, 2]
    # quadratic-form multipliers of beta, 6 x n_refl
    self._hh = numpy.array([h * h, k * k, l * l, h * k, h * l, k * l])
    # symmetry: beta_full (phasertng convention) = M . beta_ind
    c = sgtbx.tensor_rank_2_constraints(
      space_group=miller_set.space_group(), reciprocal_space=True)
    n_ind = c.n_independent_params()
    cols = []
    for i in range(n_ind):
      e = [0.] * n_ind
      e[i] = 1.
      cols.append(list(c.all_params(independent_params=e)))
    self._M = _TO_PHASERTNG.dot(numpy.array(cols).T)  # 6 x n_ind
    self._n_ind = n_ind
    self._M_pinv = numpy.linalg.pinv(self._M)
    # isotropic part removal: P = I - iso_beta (outer) iso_coeffs
    rp = uc.reciprocal_parameters()
    ap, bp, cp = rp[0], rp[1], rp[2]
    cos_r = [math.cos(math.radians(v)) for v in rp[3:6]]  # alpha*, beta*, gamma*
    iso_beta = numpy.array([ap * ap / 4, bp * bp / 4, cp * cp / 4,
      ap * bp * cos_r[2] / 2, ap * cp * cos_r[1] / 2, bp * cp * cos_r[0] / 2])
    p = uc.parameters()
    cos_d = [math.cos(math.radians(v)) for v in p[3:6]]
    iso_coeffs = (4. / 3.) * numpy.array([p[0] ** 2, p[1] ** 2, p[2] ** 2,
      p[0] * p[1] * cos_d[2], p[0] * p[2] * cos_d[1], p[1] * p[2] * cos_d[0]])
    self._P = numpy.eye(6) - numpy.outer(iso_beta, iso_coeffs)
    # sphericity restraint sigmas (RefineANISO::set_sphericity_restraint)
    sigma = 10. * (miller_set.d_min() / 2.5) ** 2
    self._sigma_sph = 0.25 * sigma * numpy.array([ap * ap, bp * bp, cp * cp,
      ap * bp, ap * cp, bp * cp])
    self._mask = (abs(self._M).sum(axis=1) > 0) & (self._sigma_sph > 0)
    self._layout()
    self._build_restraint()

  # -- parameter layout
  def _layout(self):
    self._i_bins = 1
    self._i_beta = 1 + self.n_bins
    self._i_beta1 = self._i_beta + self._n_ind
    self._n_par = self._i_beta + self._n_ind * (2 if self.resolution_dependent else 1)

  def n_parameters(self):
    return self._n_par

  def n_beta_ind(self):
    """ number of independent anisotropy parameters (of beta, and of beta1) """
    return self._n_ind

  def independent_from_beta(self, beta):
    """ independent parameters whose beta is closest to the 6-vector beta
    (phasertng order and convention) """
    return self._M_pinv.dot(numpy.asarray(beta, dtype=float))

  def beta(self, x):
    """ the 6 elements of beta, in phasertng's order and convention """
    return self._M.dot(x[self._i_beta:self._i_beta + self._n_ind])

  def beta1(self, x):
    """ the 6 elements of beta1 (zeros without resolution dependence) """
    if not self.resolution_dependent: return numpy.zeros(6)
    return self._M.dot(x[self._i_beta1:self._i_beta1 + self._n_ind])

  def start(self, i_obs, selection=None):
    """ starting parameters from the intensities (expected intensity is
    epsilon times SigmaN): the scale from the mean of I/(epsilon BEST),
    the bin corrections from the bin means relative to it (limited to a
    factor of e^3), no anisotropy. selection (numpy or flex bool, optional):
    the reflections to use """
    i = _as_numpy(i_obs)
    eps = self._eps
    ratio = i / (eps * numpy.exp(self._ln_best))
    bin_index = self._bin_index
    if selection is not None:
      sel = _as_numpy(selection, dtype=bool)
      ratio, bin_index = ratio[sel], bin_index[sel]
    scale = max(ratio.mean(), 1e-10)
    x = numpy.zeros(self._n_par)
    x[0] = math.log(scale)
    sums = numpy.bincount(bin_index, weights=ratio, minlength=self.n_bins)
    counts = numpy.bincount(bin_index, minlength=self.n_bins)
    for b in range(self.n_bins):
      if counts[b] > 0:
        m = max(sums[b] / counts[b], 1e-10)
        x[self._i_bins + b] = min(3., max(-3., math.log(m / scale)))
    return x

  # -- model
  def ln_sigman(self, x):
    """ ln SigmaN (without epsilon) per reflection, flex.double """
    x = numpy.asarray(x, dtype=float)
    v = x[0] + x[self._i_bins:self._i_bins + self.n_bins][self._bin_index] \
      + self._ln_best - self.beta(x).dot(self._hh)
    if self.resolution_dependent:
      v = v - self._x_res * (self._P.dot(self.beta1(x)).dot(self._hh))
    return flex.double(v)

  def ln_best(self):
    """ ln of the BEST curve per reflection, numpy """
    return self._ln_best

  def design(self):
    """ d ln SigmaN / dx in compact form: group_args(bin_index (numpy int
    per reflection; d/d ln_b_bin is 1 for the reflection's bin, at
    parameter 1 + bin_index), index (numpy int: the other parameters),
    rows (numpy len(index) x n_refl: their derivatives)) """
    index = [0] + list(range(self._i_beta, self._n_par))
    rows = [numpy.ones((1, self._bin_index.size)), -self._M.T.dot(self._hh)]
    if self.resolution_dependent:
      rows.append(-(self._M.T.dot(self._P.T.dot(self._hh)) * self._x_res))
    return group_args(bin_index=self._bin_index,
      index=numpy.array(index, dtype=int), rows=numpy.vstack(rows))

  def beta1_index(self):
    """ index of the first beta1 parameter (None without resolution
    dependence) """
    return self._i_beta1 if self.resolution_dependent else None

  def beta1_effective_basis(self):
    """ orthonormal basis (numpy n_beta_ind x r) of the independent beta1
    parameters that change SigmaN: P removes the isotropic part, so r is
    n_beta_ind - 1 (phasertng RefineANISO::n_resolution_dependent_parameters) """
    a = self._P.dot(self._M)
    u, s, vt = numpy.linalg.svd(a)
    r = int((s > 1e-9 * s.max()).sum()) if s.max() > 0 else 0
    return vt[:r].T

  def d_ln_sigman_dx(self, x):
    """ d ln SigmaN / dx, numpy array of shape (n_parameters, n_refl) """
    n_refl = self._bin_index.size
    d = self.design()
    g = numpy.zeros((self._n_par, n_refl))
    g[self._i_bins + self._bin_index, numpy.arange(n_refl)] = 1.
    g[d.index] = d.rows
    return g

  def anisotropy_term(self, x):
    """ exp(-h.beta_aniso.h) per reflection, the isotropic part of beta (and
    of beta1, times x) removed: nacelle's ANISOBETA column. flex.double """
    x = numpy.asarray(x, dtype=float)
    e = self._P.dot(self.beta(x)).dot(self._hh)
    if self.resolution_dependent:
      e = e + self._x_res * (self._P.dot(self.beta1(x)).dot(self._hh))
    return flex.double(numpy.exp(-e))

  # -- restraints
  def _build_restraint(self):
    """ the restraints are quadratic in x: value = x.H.x/2 """
    n = self._n_par
    h = numpy.zeros((n, n))
    w = numpy.where(self._mask, 1. / self._sigma_sph ** 2, 0.)
    e0 = numpy.zeros((6, n)); e0[:, self._i_beta:self._i_beta + self._n_ind] = self._M
    if self.resolution_dependent:
      e1 = numpy.zeros((6, n))
      e1[:, self._i_beta1:self._i_beta1 + self._n_ind] = self._M
      # anisotropy at the resolution limit, beta + P beta1, and beta1 itself
      j = self._P.dot(e0 + e1)
      h += j.T.dot(w[:, None] * j)
      j1 = self._P.dot(e1)
      h += j1.T.dot(w[:, None] * j1) / self.resolution_dependent_restraint ** 2
    else:
      j = self._P.dot(e0)
      h += j.T.dot(w[:, None] * j)
    for b in range(self.n_bins):
      i = self._i_bins + b
      h[i, i] += 1. / self._sigma_bin[b] ** 2
    self._H = h

  def restraint(self, x):
    """ (value, gradient, hessian) of the restraints (sphericity of the
    anisotropic part of beta and beta1, ln b_bin to zero), numpy """
    x = numpy.asarray(x, dtype=float)
    hx = self._H.dot(x)
    return 0.5 * x.dot(hx), hx, self._H.copy()

# ---------------------------------------------------------------------------
# Joint data-only fit of SigmaN and the tNCS coherence

normalisation_params = iotbx.phil.parse("""\
  gamma_max = 1.05
    .type = float
    .short_caption = Upper limit of the tNCS coherence gamma
    .help = "Upper bound of gamma, the low-resolution limit of the coherence " \\
            "of the tNCS-related copies in the data relative to the " \\
            "consensus of the model's copies. Values above 1 (a stronger " \\
            "modulation than identical copies give) are allowed up to this " \\
            "limit and reported with a warning."
  outlier_probability = 1e-6
    .type = float
    .short_caption = Wilson outlier probability
    .help = "Reflections whose Wilson tail probability (including the " \\
            "measurement error) is below this are excluded from the fit."
  resolution_dependent_restraint = 10.0
    .type = float
    .short_caption = Restraint factor of the resolution-dependent anisotropy
    .help = "Factor multiplying the sigma of the sphericity restraint on " \\
            "the resolution-dependent anisotropy term (beta1), relative " \\
            "to that of the anisotropy at the resolution limit; as " \\
            "phasertng macaniso.resolution_dependent_restraint."
  resolution_dependent_anisotropy = off on *auto
    .type = choice
    .short_caption = Resolution-dependent anisotropy
    .help = "Anisotropy that changes linearly with d*^2 (phasertng " \\
            "macaniso.resolution_dependent). auto: refined after the fit " \\
            "without it and kept if the Wilson log-likelihood gain exceeds " \\
            "half the 99.9% point of chi-squared for its parameters, as " \\
            "phasertng runANISO."
  bins
    .short_caption = Resolution bins of the SigmaN correction
    .help = "Bins of equal width in d*^2, as phasertng Bin::setup: the " \\
            "number of reflections divided by width, within range, and " \\
            "at least 10 reflections per bin."
  {
    width = 500
      .type = int
      .short_caption = Reflections per bin
      .help = "Target number of reflections per bin."
    range = 6 50
      .type = ints(size=2)
      .short_caption = Smallest and largest number of bins
      .help = "Smallest and largest number of bins."
  }
  hierarchical_coherence = False
    .type = bool
    .short_caption = Coherence per tNCS group
    .help = "Refine a coherence gamma_g = gamma exp(eta_g) per tNCS group, " \\
            "eta_g restrained to zero with sigma_eta, instead of one shared " \\
            "gamma."
  sigma_eta = 0.1
    .type = float
    .short_caption = Restraint sigma of ln(gamma_g/gamma)
    .help = "Standard deviation of the restraint of eta_g = ln(gamma_g/" \\
            "gamma) to zero (hierarchical_coherence)."
  max_iterations = 50
    .type = int
    .short_caption = Newton iterations per fit
    .help = "Largest number of Newton iterations per fit (each outlier " \\
            "round, and each anisotropy model, is a fit)."
""")

# phasertng run/runANISO.cc lines 146-150 (resolution_dependent == "auto"):
# chisq_999, the 99.9% points of chi-squared for 1..5 parameters; the
# threshold of the Wilson log-likelihood gain (without restraints, on the
# reflections selected by both fits) is half of it, for min(npar1, 5)
# parameters, npar1 = n_beta_ind - 1 (RefineANISO.cc line 1054)
CHISQ_999 = (10.828, 13.816, 16.266, 18.467, 20.515)

MAX_OUTLIER_ROUNDS = 3

GAMMA_WARNING = ("tNCS coherence gamma = %.3f exceeds 1: the data are more "
  "strongly modulated than identical copies would give")

def beta_as_b_cart(unit_cell, beta):
  """ the anisotropic B (Cartesian, A^2) of beta in phasertng's convention,
  exp(-(b11 h^2 + ... + b12 hk + ...)) = exp(-h.B*.h/4) """
  b = numpy.asarray(beta, dtype=float) / numpy.array([1, 1, 1, 2, 2, 2.])
  u_star = tuple(b / (2 * math.pi ** 2))
  return adptbx.u_as_b(adptbx.u_star_as_u_cart(unit_cell, u_star))

def delta_b_aniso(unit_cell, beta):
  """ largest minus smallest eigenvalue of beta_as_b_cart (A^2) """
  e = adptbx.eigenvalues(beta_as_b_cart(unit_cell, beta))
  return max(e) - min(e)

def _wilson_outliers(i, sig, sigma, centric, p_min):
  """ True where the Wilson tail probability of the observation
  (xray_ext.wilson_intensity_tail_probability; i, sig, sigma numpy with
  sig, sigma > 0) is below p_min. The tail probability (up to ~200
  quadratures per centric reflection) is evaluated only where a lower bound
  of it is below 100 p_min: with J the true and I = J + e the observed
  intensity, P(I > i) >= P(J > i) P(e >= 0) for i >= 0 and
  P(I < i) >= P(J < d) P(e < i - d) for i < 0 and any d > 0. Returns
  (numpy bool, number of reflections evaluated) """
  from scipy.special import erf, erfc, ndtr
  e = i / sigma
  s = sig / sigma
  centric = numpy.asarray(centric, dtype=bool)
  bound = numpy.zeros(e.size)
  pos = e >= 0
  pa, pc = pos & ~centric, pos & centric
  bound[pa] = 0.5 * numpy.exp(-e[pa])
  bound[pc] = 0.5 * erfc(numpy.sqrt(e[pc] / 2))
  neg = ~pos
  en, sn, cn = e[neg], s[neg], centric[neg]
  bn = numpy.zeros(en.size)
  for k in (0.1, 0.3, 1., 3., 10.):
    d = k * sn
    p_j = numpy.where(cn, erf(numpy.sqrt(d / 2)), -numpy.expm1(-d))
    bn = numpy.maximum(bn, p_j * ndtr((en - d) / sn))
  bound[neg] = bn
  isel = numpy.flatnonzero(bound < 100 * p_min)
  flags = numpy.zeros(e.size, dtype=bool)
  if isel.size:
    p = xray_ext.wilson_intensity_tail_probability(
      i_obs=flex.double(i[isel]), sig_i_obs=flex.double(sig[isel]),
      sigma=flex.double(sigma[isel]),
      centric_flags=flex.bool(centric[isel])).as_numpy_array()
    flags[isel] = p < p_min
  return flags, isel.size

def _fit_data(i_obs, tncs_a_groups, params):
  """ the data of the fit as numpy arrays, with the resolution bins and the
  tNCS pair terms (numpy n_groups x n_refl, or None) """
  assert i_obs.sigmas() is not None
  assert not i_obs.anomalous_flag()
  i = _as_numpy(i_obs.data())
  sig = _as_numpy(i_obs.sigmas())
  d_star_sq = i_obs.d_star_sq().data()
  dss = _as_numpy(d_star_sq)
  present = numpy.isfinite(i) & numpy.isfinite(sig) & (sig > 0)
  # I == 0 exactly is missing (as in nacelle)
  zero_i = present & (i == 0)
  valid = present & ~zero_i
  bins = resolution_bins(d_star_sq, width=params.bins.width,
    n_min=params.bins.range[0], n_max=params.bins.range[1])
  a = None
  if tncs_a_groups:
    a = numpy.vstack([_as_numpy(g) for g in tncs_a_groups])
    assert a.shape[1] == i.size
  return group_args(i=i, sig=sig, valid=valid, zero_i=zero_i, bins=bins, a=a,
    eps=_as_numpy(i_obs.epsilons().data()),
    centric=_as_numpy(i_obs.centric_flags().data(), dtype=bool),
    c=2 * math.pi ** 2 * dss)

class _target(object):
  """ minus the Wilson log-likelihood of the selected intensities plus the
  restraints, as a function of p = [q, gamma, w, eta_g...]: q are the
  parameters of the sigman_model with its beta1 reduced to the components
  that change SigmaN (beta1_effective_basis), w = v^2 (so that the
  coherence is smooth at v = 0), and eta_g = ln(gamma_g/gamma) with
  hierarchical coherence. The expected intensity is
  Sigma = epsilon SigmaN TEPS, TEPS = max(1 + gamma exp(-c w) Z, TEPS_FLOOR),
  c = 2 pi^2 s^2, Z = sum_g exp(eta_g) A_g. Derivatives are by the chain
  rule over the per-reflection derivatives in ln Sigma """

  def __init__(self, data, model, selection, hierarchical=False,
      sigma_eta=0.1):
    self.data, self.model = data, model
    sel = numpy.asarray(selection, dtype=bool) & data.valid
    self.selection = sel
    isel = numpy.flatnonzero(sel)
    d = model.design()
    n_x = model.n_parameters()
    # x = E q: identity except for beta1 = Q y
    i1 = model.beta1_index()
    if i1 is None:
      e = numpy.eye(n_x)
    else:
      nb1 = model.n_beta_ind()
      basis = model.beta1_effective_basis()
      e = numpy.zeros((n_x, i1 + basis.shape[1]))
      e[:i1, :i1] = numpy.eye(i1)
      e[i1:i1 + nb1, i1:] = basis
    self.expansion = e
    self.n_q = e.shape[1]
    self.n_bins = model.n_bins
    # dense rows in q (all parameters but the bins)
    rows_x = numpy.zeros((n_x, d.bin_index.size))
    rows_x[d.index] = d.rows
    dense_q = [0] + list(range(1 + self.n_bins, self.n_q))
    self.dense_index = numpy.array(dense_q, dtype=int)
    self.rows = e.T[self.dense_index].dot(rows_x)[:, isel]
    del rows_x
    self.bin = d.bin_index[isel]
    self.h_q = e.T.dot(model.restraint(numpy.zeros(n_x))[2]).dot(e)
    self.ln_base = (numpy.log(data.eps) + model.ln_best())[isel]
    self.i_np, self.sig_np = data.i[isel], data.sig[isel]
    self.cen_np = data.centric[isel]
    self.f_i, self.f_sig = flex.double(self.i_np), flex.double(self.sig_np)
    self.f_cen = flex.bool(self.cen_np)
    self.fisher_k = numpy.where(self.cen_np, 2., 1.)
    self.isel = isel
    # tNCS
    self.tncs = data.a is not None
    self.hierarchical = bool(hierarchical) and self.tncs
    self.sigma_eta = sigma_eta
    self.i_tncs = self.n_q
    self.n_groups = data.a.shape[0] if self.tncs else 0
    n_t = 0
    if self.tncs:
      n_t = 2 + (self.n_groups if self.hierarchical else 0)
      self.a = data.a[:, isel]
      self.a_sum = self.a.sum(axis=0)
      self.c = data.c[isel]
    self.n_t = n_t
    self.n_p = self.n_q + n_t

  def n_selected(self):
    return self.isel.size

  # -- parameters
  def start(self, i_obs):
    """ sigman_model.start on the selected reflections, gamma = 1, v = 0,
    eta = 0 """
    p = numpy.zeros(self.n_p)
    x = self.model.start(i_obs.data(), selection=self.selection)
    p[:self.n_q] = self.reduce(x)
    if self.tncs: p[self.i_tncs] = 1.
    return p

  def reduce(self, x):
    """ q (least squares) from the sigman_model x """
    return numpy.linalg.lstsq(self.expansion, x, rcond=None)[0]

  def expand(self, p):
    """ the sigman_model x of p """
    return self.expansion.dot(p[:self.n_q])

  def bounds(self, gamma_max):
    lower = numpy.full(self.n_p, -numpy.inf)
    upper = numpy.full(self.n_p, numpy.inf)
    if self.tncs:
      lower[self.i_tncs] = 0.
      upper[self.i_tncs] = gamma_max
      lower[self.i_tncs + 1] = 0.
    return lower, upper

  def gamma_groups(self, p):
    if not self.hierarchical: return None
    return [float(v) for v in p[self.i_tncs] * numpy.exp(p[self.i_tncs + 2:])]

  # -- tNCS modulation
  def _teps_raw(self, p, a, a_sum, c):
    gamma, w = p[self.i_tncs], p[self.i_tncs + 1]
    g = numpy.exp(-c * w)
    if self.hierarchical:
      eg = numpy.exp(p[self.i_tncs + 2:])
      z = eg.dot(a)
    else:
      eg, z = None, a_sum
    return 1 + gamma * g * z, g, z, eg

  def teps_all(self, p, floor=None):
    """ TEPS for all reflections (numpy) and the number floored, at U_MIN
    (as reported) or at floor """
    d = self.data
    if floor is None: floor = U_MIN
    if not self.tncs: return numpy.ones(d.i.size), 0
    t = self._teps_raw(p, d.a, d.a.sum(axis=0), d.c)[0]
    low = t < floor
    return numpy.where(low, floor, t), int(low.sum())

  def ln_sigman_all(self, p):
    return _as_numpy(self.model.ln_sigman(self.expand(p)))

  def sigma_all(self, p):
    """ epsilon SigmaN TEPS for all reflections (numpy), with TEPS floored
    at TEPS_FLOOR as in the target """
    return self.data.eps * numpy.exp(self.ln_sigman_all(p)) * \
      self.teps_all(p, TEPS_FLOOR)[0]

  # -- target
  def _restraint(self, p):
    q = p[:self.n_q]
    hq = self.h_q.dot(q)
    v = 0.5 * q.dot(hq)
    g = numpy.zeros(self.n_p)
    h = numpy.zeros((self.n_p, self.n_p))
    g[:self.n_q] = hq
    h[:self.n_q, :self.n_q] = self.h_q
    if self.hierarchical:
      eta = p[self.i_tncs + 2:]
      w = 1. / self.sigma_eta ** 2
      v += 0.5 * w * (eta * eta).sum()
      j = numpy.arange(self.i_tncs + 2, self.n_p)
      g[j] = w * eta
      h[j, j] = w
    return v, g, h

  def _ln_sigma(self, p):
    q = p[:self.n_q]
    ln_sigma = self.ln_base + q[self.dense_index].dot(self.rows) + \
      q[1 + self.bin]
    tn = None
    if self.tncs:
      t, g, z, eg = self._teps_raw(p, self.a, self.a_sum, self.c)
      low = t < TEPS_FLOOR
      tn = group_args(t=numpy.where(low, TEPS_FLOOR, t), g=g, z=z, eg=eg,
        low=low)
      ln_sigma = ln_sigma + numpy.log(tn.t)
    return ln_sigma, tn

  def _wilson(self, ln_sigma):
    if not numpy.isfinite(ln_sigma).all() or abs(ln_sigma).max() > 700:
      return None
    sigma = numpy.exp(ln_sigma)
    w = xray_ext.wilson_intensity(i_obs=self.f_i, sig_i_obs=self.f_sig,
      sigma=flex.double(sigma), centric_flags=self.f_cen)
    ll = w.log_likelihood.as_numpy_array()
    llsum = ll.sum()
    if not math.isfinite(llsum): return None
    return group_args(w=w, sigma=sigma, ll=llsum)

  def log_likelihood(self, p):
    """ the Wilson log-likelihood (no restraints), -inf if not finite """
    r = self._wilson(self._ln_sigma(p)[0])
    return -numpy.inf if r is None else r.ll

  def value(self, p):
    """ the target (inf where Sigma is out of range) """
    r = self._wilson(self._ln_sigma(p)[0])
    if r is None: return numpy.inf
    return -r.ll + self._restraint(p)[0]

  def evaluate(self, p, hessian="exact"):
    """ (target, gradient, hessian or None); hessian: "exact", "fisher"
    (sum over the reflections of w g g^T, with the Fisher information
    w = 1/(k + (sig_i/Sigma)^2) in ln Sigma, k = 1 acentric and 2 centric,
    plus the restraints) or None """
    ln_sigma, tn = self._ln_sigma(p)
    r = self._wilson(ln_sigma)
    if r is None: return numpy.inf, None, None
    rv, rg, rh = self._restraint(p)
    f = -r.ll + rv
    d1 = r.w.d_ll_d_ln_sigma.as_numpy_array()
    nb = self.n_bins
    rows = self.rows
    kidx = self.dense_index
    if self.tncs:
      # d ln TEPS / d (gamma, w, eta_g); zero where TEPS is floored
      gz = numpy.where(tn.low, 0., tn.g * tn.z / tn.t)
      gamma = p[self.i_tncs]
      rt = [gz, -self.c * gamma * gz]
      if self.hierarchical:
        ga = numpy.where(tn.low, 0., tn.g / tn.t) * tn.eg[:, None] * self.a
        rt.extend(list(gamma * ga))
      rt = numpy.vstack(rt)
      rows = numpy.vstack([rows, rt])
      kidx = numpy.concatenate([kidx, numpy.arange(self.i_tncs, self.n_p)])
    g_ll = numpy.zeros(self.n_p)
    g_ll[kidx] = rows.dot(d1)
    g_ll[1:1 + nb] = numpy.bincount(self.bin, weights=d1, minlength=nb)
    grad = -g_ll + rg
    if hessian is None: return f, grad, None
    if hessian == "exact":
      wt = r.w.d2_ll_d_ln_sigma2.as_numpy_array()
    else:
      assert hessian == "fisher"
      s = self.sig_np / r.sigma
      wt = -1. / (self.fisher_k + s * s)
    h_ll = numpy.zeros((self.n_p, self.n_p))
    h_ll[numpy.ix_(kidx, kidx)] = (rows * wt).dot(rows.T)
    ib = numpy.arange(1, 1 + nb)
    h_ll[ib, ib] = numpy.bincount(self.bin, weights=wt, minlength=nb)
    for j, k in enumerate(kidx):
      col = numpy.bincount(self.bin, weights=wt * rows[j], minlength=nb)
      h_ll[ib, k] = col
      h_ll[k, ib] = col
    if hessian == "exact" and self.tncs:
      # the second derivatives of ln TEPS, sum of d1 (d2 T/T - r r^T)
      it = numpy.arange(self.i_tncs, self.n_p)
      n_t = self.n_t
      d2t = numpy.zeros((n_t, n_t))
      c = self.c
      d2t[0, 1] = d2t[1, 0] = (d1 * (-c * rt[0])).sum()
      d2t[1, 1] = (d1 * (-c * rt[1])).sum()
      if self.hierarchical:
        for gi in range(self.n_groups):
          j = 2 + gi
          d2t[0, j] = d2t[j, 0] = (d1 * ga[gi]).sum()
          d2t[1, j] = d2t[j, 1] = (d1 * (-c * rt[j])).sum()
          d2t[j, j] = (d1 * rt[j]).sum()
      d2t -= (rt * d1).dot(rt.T)
      h_ll[numpy.ix_(it, it)] += d2t
    return f, grad, -h_ll + rh

def _newton_step(h, g, free, ridge=False):
  """ -H^-1 g on the free parameters, scaled to unit diagonal; None if H is
  not positive definite there (and ridge is False). With ridge, a multiple
  of the identity is added until it is """
  step = numpy.zeros(g.size)
  j = numpy.flatnonzero(free)
  if j.size == 0: return step
  hs = h[numpy.ix_(j, j)]
  d = numpy.sqrt(abs(numpy.diag(hs)))
  d[d == 0] = 1.
  hs = hs / numpy.outer(d, d)
  gs = g[j] / d
  for lam in ((0.,) + ((1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1., 1e2) if ridge
      else ())):
    try:
      l = numpy.linalg.cholesky(hs + lam * numpy.eye(j.size))
    except numpy.linalg.LinAlgError:
      continue
    y = numpy.linalg.solve(l, gs)
    step[j] = -numpy.linalg.solve(l.T, y) / d
    return step
  return None

def _minimise(t, p, lower, upper, max_iterations, tolerance=1e-4):
  """ Newton steps with the exact Hessian (the Fisher form where that is not
  positive definite for the target, i.e. the Hessian of the likelihood not
  negative definite), bounds by an active set and projection, and a
  backtracking line search on the target. Returns group_args(p, f,
  converged, n_iterations, n_fisher) """
  p = numpy.clip(p, lower, upper)
  f, g, h = t.evaluate(p)
  assert math.isfinite(f)
  converged = False
  n_fisher = 0
  it = 0
  for it in range(1, max_iterations + 1):
    free = ~(((p <= lower) & (g > 0)) | ((p >= upper) & (g < 0)))
    step = _newton_step(h, g, free)
    if step is None:
      n_fisher += 1
      step = _newton_step(t.evaluate(p, hessian="fisher")[2], g, free,
        ridge=True)
    decrement = -g.dot(step)
    if decrement < tolerance:
      converged = True
      break
    alpha, accepted = 1., False
    for k in range(40):
      pn = numpy.clip(p + alpha * step, lower, upper)
      fn = t.value(pn)
      if fn <= f + 1e-4 * g.dot(pn - p):
        accepted = True
        break
      alpha *= 0.5
    if not accepted:
      converged = decrement < 100 * tolerance
      break
    p = pn
    f, g, h = t.evaluate(p)
  return group_args(p=p, f=f, converged=converged, n_iterations=it,
    n_fisher=n_fisher)

def _fit_model(i_obs, data, params, resolution_dependent, selection=None,
    start=None, start_selection=None):
  """ one sigman_model (with or without beta1) and the coherence, with
  outlier rounds unless selection is given. start: p of a fit (its q
  extended with zeros for beta1) """
  model = sigman_model(i_obs, data.bins,
    resolution_dependent=resolution_dependent,
    resolution_dependent_restraint=params.resolution_dependent_restraint)
  fixed = selection is not None
  if fixed:
    sel = _as_numpy(selection, dtype=bool) & data.valid
  elif start_selection is not None:
    sel = start_selection & data.valid
  else:
    sel = data.valid.copy()
  def target(sel):
    return _target(data, model, sel, hierarchical=params.hierarchical_coherence,
      sigma_eta=params.sigma_eta)
  t = target(sel)
  if start is None:
    p = t.start(i_obs)
  else:
    p = numpy.zeros(t.n_p)
    n_q0 = start.n_q
    p[:n_q0] = start.p[:n_q0]
    p[t.n_q:] = start.p[n_q0:]
  lower, upper = t.bounds(params.gamma_max)
  n_rounds, n_iterations, n_fisher = 0, 0, 0
  converged = True
  iv = numpy.flatnonzero(data.valid)
  while True:
    m = _minimise(t, p, lower, upper, params.max_iterations)
    p = m.p
    n_iterations += m.n_iterations
    n_fisher += m.n_fisher
    converged = converged and m.converged
    sigma = t.sigma_all(p)
    out = _wilson_outliers(data.i[iv], data.sig[iv], sigma[iv],
      data.centric[iv], params.outlier_probability)[0]
    new_sel = numpy.zeros(data.valid.size, dtype=bool)
    new_sel[iv] = ~out
    n_change = int((new_sel != sel).sum())
    if fixed or n_change == 0 or n_rounds == MAX_OUTLIER_ROUNDS: break
    n_rounds += 1
    sel = new_sel
    t = target(sel)
  return group_args(model=model, target=t, p=p, n_q=t.n_q, f=m.f,
    selection=sel, n_outlier_rounds=n_rounds, n_would_change=n_change,
    n_iterations=n_iterations, n_fisher=n_fisher, converged=converged,
    log_likelihood=t.log_likelihood(p))

def fit(i_obs, tncs_a_groups=None, params=None, log=None, selection=None):
  """ SigmaN (sigman_model: BEST, bin corrections, anisotropy) and, with
  tNCS, the coherence Gamma(s) = gamma exp(-2 pi^2 v^2 s^2) shared by the
  groups, by the exact Wilson likelihood of the intensities with their
  measurement errors. Newton steps (exact Hessian, the Fisher form where
  that is not definite) with a backtracking line search, 0 <= gamma <=
  gamma_max and v >= 0, starting at gamma = 1, v = 0. Reflections with
  SIGI <= 0, and with I == 0 exactly (as nacelle), are excluded; outliers (Wilson tail probability below
  outlier_probability) are excluded and the fit repeated, at most 3 times,
  unless selection (flex.bool) is given, which is used as is (n_would_change
  then counts the reflections whose outlier status would differ).

  i_obs: non-anomalous intensities with sigmas, in the asu.
  tncs_a_groups: None, or a list of flex.double, the consensus pair term
    A_g of each tNCS group on the indices of i_obs.
  Returns group_args with resn (sqrt(epsilon SigmaN)), teps
  (max(1 + Gamma sum_g A_g, U_MIN), or 1; the fit itself floors TEPS at
  TEPS_FLOOR), anisotropy (miller.arrays),
  selection (flex.bool, True = used), coherence (coherence_curve or None),
  gamma_groups (hierarchical_coherence: gamma_g per group, else None), x
  (sigman_model parameters), model, log_likelihood (Wilson, selected, no
  restraints), n_outliers, n_zero_sigma, n_zero_intensity, n_teps_floored, n_would_change,
  n_outlier_rounds, resolution_dependent (bool) with its gain and
  threshold, delta_b_aniso (A^2, of beta), converged, warnings. """
  if params is None: params = normalisation_params.extract()
  if log is None: log = null_out()
  if not tncs_a_groups: tncs_a_groups = None
  data = _fit_data(i_obs, tncs_a_groups, params)
  n_zero_intensity = int(data.zero_i.sum())
  n_zero_sigma = int((~data.valid).sum()) - n_zero_intensity
  mode = params.resolution_dependent_anisotropy
  r = _fit_model(i_obs, data, params, resolution_dependent=(mode == "on"),
    selection=selection)
  gain = threshold = None
  if mode == "auto":
    n1 = r.model.n_beta_ind() - 1
    if n1 > 0:
      r1 = _fit_model(i_obs, data, params, resolution_dependent=True,
        selection=selection, start=r, start_selection=r.selection)
      # both likelihoods on the reflections selected by both (each model's
      # outlier rejection would favour that model), as runANISO
      both = r.selection & r1.selection
      def ll(rr):
        return _target(data, rr.model, both,
          hierarchical=params.hierarchical_coherence,
          sigma_eta=params.sigma_eta).log_likelihood(rr.p)
      gain = ll(r1) - ll(r)
      threshold = CHISQ_999[min(n1, 5) - 1] / 2.
      print("Resolution-dependent anisotropy: Wilson log-likelihood gain "
        "%.2f for %d parameter(s), threshold %.2f: %s" % (gain, n1,
        threshold, "used" if gain > threshold else "not used"), file=log)
      if gain > threshold: r = r1
  t, p, model = r.target, r.p, r.model
  x = t.expand(p)
  sel = r.selection
  eps = data.eps
  resn = numpy.sqrt(eps * numpy.exp(t.ln_sigman_all(p)))
  teps, n_floored = t.teps_all(p)
  coherence = gamma_groups = None
  warnings = []
  if t.tncs:
    gamma, w = p[t.i_tncs], p[t.i_tncs + 1]
    coherence = coherence_curve(gamma, math.sqrt(max(w, 0.)))
    gamma_groups = t.gamma_groups(p)
    g_max = max([gamma] + (gamma_groups or []))
    if g_max > 1: warnings.append(GAMMA_WARNING % g_max)
  if not r.converged:
    warnings.append("normalisation fit stopped at max_iterations = %d "
      "before convergence" % params.max_iterations)
  def array(v):
    return i_obs.array(data=flex.double(v))
  n_outliers = int((data.valid & ~sel).sum())
  result = group_args(
    resn=array(resn),
    teps=array(teps),
    anisotropy=i_obs.array(data=model.anisotropy_term(x)),
    selection=flex.bool(sel),
    coherence=coherence,
    gamma_groups=gamma_groups,
    x=x,
    model=model,
    log_likelihood=r.log_likelihood,
    n_outliers=n_outliers,
    n_zero_sigma=n_zero_sigma,
    n_zero_intensity=n_zero_intensity,
    n_teps_floored=n_floored,
    n_would_change=r.n_would_change,
    n_outlier_rounds=r.n_outlier_rounds,
    n_iterations=r.n_iterations,
    resolution_dependent=model.resolution_dependent,
    resolution_dependent_gain=gain,
    resolution_dependent_threshold=threshold,
    delta_b_aniso=delta_b_aniso(i_obs.unit_cell(), model.beta(x)),
    converged=r.converged,
    warnings=warnings)
  print("LLGI data normalisation: %d reflections used, %d outliers, %d "
    "with SIGI <= 0, %d with I = 0" % (sel.sum(), n_outliers,
    n_zero_sigma, n_zero_intensity), file=log)
  print("  Wilson log-likelihood %.2f, %d Newton iterations, %d outlier "
    "round(s), anisotropy delta B %.2f A^2%s" % (r.log_likelihood,
    r.n_iterations, r.n_outlier_rounds, result.delta_b_aniso,
    " (resolution-dependent)" if model.resolution_dependent else ""),
    file=log)
  if coherence is not None:
    print("  tNCS %s%s; %d reflections with TEPS at %.2f" % (
      coherence.summary(), "" if gamma_groups is None else
      ", per group " + " ".join(["%.3f" % v for v in gamma_groups]),
      n_floored, U_MIN), file=log)
  for w in warnings: print("  Warning: " + w, file=log)
  return result

# ---------------------------------------------------------------------------
# Patterson test for translational NCS

def _origin_distance_of(unit_cell, patterson_ops, translations, site):
  """ (distance, site) of the position closest to the origin among the
  symmetry equivalents of site (the operators of the Patterson group), with
  the lattice translations of the space group (those of its centring) and
  of the unit cell added """
  best = None
  for op in patterson_ops:
    s0 = op * site
    for lt in translations:
      for i in (-1, 0, 1):
        for j in (-1, 0, 1):
          for k in (-1, 0, 1):
            s = (s0[0] + lt[0] + i, s0[1] + lt[1] + j, s0[2] + lt[2] + k)
            d = unit_cell.length(s)
            if(best is None or d < best[0]): best = (d, s)
  return best

def patterson_tncs_check(i_obs, selection=None, resolution=(5.0, 15.0),
      percent=16., origin_distance=15., minimum_resolution=8.,
      minimum_reflections=50):
  """ Whether the intensities have translational NCS, by the peak test of
  phasertng's tncso (run/runTNCSO.cc, io/PattersonFunction.cc): the
  Patterson function of the positive intensities (in the symmetry of the
  Patterson group of the space group with the screw and glide translations
  removed and the centring kept, build_derived_group(False, False) as
  phasertng) of
  the selected reflections with d between max(d_min, resolution[0]) and
  min(d_max, resolution[1]) (all the selected ones if none are), and a peak
  search of it (maptbx.peak_search, level 3, interpolated). tNCS is
  indicated if a peak has a height of at least percent (% of the origin
  peak) at a distance of at least origin_distance (A) from the origin and
  from the lattice points equivalent to the origin (symmetry equivalents
  considered, with the cell and centring translations).
  No decision (tncs = None) if there are fewer than minimum_reflections
  reflections, the Patterson has a resolution worse than minimum_resolution
  + 0.2 (as phasertng, for the sampling of Miller indices), or a cell edge
  is shorter than origin_distance.

  i_obs: intensities (anomalous ones are averaged); selection: flex.bool
  on its indices, or None.
  Returns group_args(tncs (True, False or None), peaks (list of
  (height_percent, distance_A, site_frac) of the non-origin peaks above
  the threshold, highest first, the site the equivalent closest to the
  origin), d_min, d_max (of the Patterson coefficients; None if there are
  none), top ((height_percent, distance_A, site_frac) of the highest
  non-origin peak beyond origin_distance, or None), reason (str)). """
  from cctbx import crystal, miller, maptbx
  if(i_obs.anomalous_flag()):
    i_obs = i_obs.average_bijvoet_mates()
  if(selection is not None):
    i_obs = i_obs.select(selection)
  i_obs = i_obs.map_to_asu()
  def result(tncs, reason, peaks=(), d_min=None, d_max=None, top=None):
    return group_args(tncs=tncs, peaks=list(peaks), d_min=d_min,
      d_max=d_max, reason=reason, top=top)
  if(i_obs.size() == 0):
    return result(None, "no reflections")
  unit_cell = i_obs.unit_cell()
  d = i_obs.d_spacings().data()
  data = i_obs.data()
  hires = max(flex.min(d), resolution[0])
  lores = min(flex.max(d), resolution[1])
  use = (data > 0) & (d >= hires - 1e-6) & (d <= lores + 1e-6)
  if(use.count(True) == 0):
    use = data > 0 # as phasertng: if none, all the reflections
  n = use.count(True)
  if(n < minimum_reflections):
    return result(None, "only %d reflections for the Patterson (at least "
      "%d needed)" % (n, minimum_reflections))
  coefficients = i_obs.select(use)
  d_max, d_min = coefficients.d_max_min()
  if(d_min > minimum_resolution + 0.2):
    return result(None, "the resolution of the Patterson, %.1f A, is worse "
      "than %.1f A" % (d_min, minimum_resolution), d_min=d_min, d_max=d_max)
  if(min(unit_cell.parameters()[:3]) < origin_distance):
    return result(None, "a cell edge is shorter than %.1f A" %
      origin_distance, d_min=d_min, d_max=d_max)
  space_group = i_obs.space_group()
  derived = space_group.build_derived_group(False, False)
  derived_symmetry = crystal.symmetry(unit_cell=unit_cell,
    space_group=derived)
  f_patt = miller.array(
    miller.set(derived_symmetry, coefficients.indices(),
      anomalous_flag=False),
    data=flex.sqrt(coefficients.data())).map_to_asu()
  patterson_map = f_patt.patterson_map(d_min=d_min,
    symmetry_flags=maptbx.use_space_group_symmetry)
  peak_list = patterson_map.peak_search(
    maptbx.peak_search_parameters(peak_search_level=3, max_peaks=0,
      interpolate=True), verify_symmetry=False)
  heights = peak_list.heights()
  if(heights.size() == 0):
    return result(False, "no peaks in the Patterson", d_min=d_min,
      d_max=d_max)
  origin_height = flex.max(heights)
  patterson_ops = list(derived.build_derived_patterson_group().all_ops())
  translations = [tuple(lt.as_double()) for lt in space_group.ltr()]
  peaks = []
  top = None
  for height, site in zip(heights, peak_list.sites()):
    height_percent = 100. * height / origin_height
    if(top is not None and height_percent < percent):
      break # the peak list is sorted by height; top is known
    closest = _origin_distance_of(unit_cell, patterson_ops, translations,
      site)
    if(closest[0] < origin_distance): continue
    peak = (height_percent, closest[0], closest[1])
    if(top is None or height_percent > top[0]): top = peak
    if(height_percent >= percent): peaks.append(peak)
  peaks.sort(key=lambda p: -p[0])
  if(len(peaks) > 0):
    reason = ("a Patterson peak of %.1f%% of the origin peak at %.1f A "
      "(at least %.1f%% beyond %.1f A needed)" % (peaks[0][0], peaks[0][1],
        percent, origin_distance))
  else:
    reason = "no Patterson peak of %.1f%% or more of the origin peak " \
      "beyond %.1f A" % (percent, origin_distance)
  return result(len(peaks) > 0, reason, peaks, d_min, d_max, top)
