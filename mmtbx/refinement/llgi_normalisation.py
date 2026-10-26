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
"""
import math
import numpy
from libtbx import group_args
from cctbx import sgtbx
from cctbx.array_family import flex

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

# phasertng RefineANISO: RESOLUTION_DEPENDENT_SIGMA_FACTOR
RESOLUTION_DEPENDENT_SIGMA_FACTOR = 1.
# phasertng RefineANISO: LOG_WILSON_SIGMA_CONSTANT, LOG_WILSON_SIGMA_FACTOR
LOG_WILSON_SIGMA_CONSTANT = 0.15
LOG_WILSON_SIGMA_FACTOR = 0.001

_BEST_TABLE_NP = numpy.array(BEST_TABLE)
assert len(BEST_TABLE) == 300

def best_curve_values(d_star_sq):
  """ BEST curve (times k) at d*^2, as phasertng's best_curve_at_ssqr::get:
  linear interpolation, truncating the table index towards zero, constant
  for an index below 0 and from 299 on. flex.double in, flex.double out """
  ssqr = numpy.array(list(d_star_sq), dtype=float)
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
  return flex.double(list(BEST_K * v))

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
  ssqr = numpy.array(list(d_star_sq), dtype=float)
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
  return group_args(index=flex.int(list(index)), n_bins=n,
    d_mid=flex.double(list(d_mid)))

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

  def __init__(self, miller_set, bins, resolution_dependent=False):
    self.n_bins = bins.n_bins
    self.resolution_dependent = bool(resolution_dependent)
    self._bin_index = numpy.array(list(bins.index), dtype=int)
    self._eps = numpy.array(list(miller_set.epsilons().data()), dtype=float)
    self._sigma_bin = numpy.array(list(bin_restraint_sigmas(bins.d_mid)))
    uc = miller_set.unit_cell()
    ssqr = numpy.array(list(miller_set.d_star_sq().data()))
    self._ln_best = numpy.log(numpy.array(list(best_curve_values(
      flex.double(list(ssqr))))))
    ssqr_max = ssqr.max()
    self._x_res = ssqr / ssqr_max
    hkl = numpy.array(list(miller_set.indices()), dtype=float)
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

  def start(self, i_obs):
    """ starting parameters from the intensities (expected intensity is
    epsilon times SigmaN): the scale from the mean of I/(epsilon BEST),
    the bin corrections from the bin means relative to it (limited to a
    factor of e^3), no anisotropy """
    i = numpy.array(list(i_obs), dtype=float)
    eps = self._eps
    ratio = i / (eps * numpy.exp(self._ln_best))
    scale = max(ratio.mean(), 1e-10)
    x = numpy.zeros(self._n_par)
    x[0] = math.log(scale)
    sums = numpy.bincount(self._bin_index, weights=ratio, minlength=self.n_bins)
    counts = numpy.bincount(self._bin_index, minlength=self.n_bins)
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
    return flex.double(list(v))

  def d_ln_sigman_dx(self, x):
    """ d ln SigmaN / dx, numpy array of shape (n_parameters, n_refl) """
    n_refl = self._bin_index.size
    g = numpy.zeros((self._n_par, n_refl))
    g[0] = 1.
    g[self._i_bins + self._bin_index, numpy.arange(n_refl)] = 1.
    g[self._i_beta:self._i_beta + self._n_ind] = -self._M.T.dot(self._hh)
    if self.resolution_dependent:
      g[self._i_beta1:self._i_beta1 + self._n_ind] = -(
        self._M.T.dot(self._P.T.dot(self._hh)) * self._x_res)
    return g

  def anisotropy_term(self, x):
    """ exp(-h.beta_aniso.h) per reflection, the isotropic part of beta (and
    of beta1, times x) removed: nacelle's ANISOBETA column. flex.double """
    x = numpy.asarray(x, dtype=float)
    e = self._P.dot(self.beta(x)).dot(self._hh)
    if self.resolution_dependent:
      e = e + self._x_res * (self._P.dot(self.beta1(x)).dot(self._hh))
    return flex.double(list(numpy.exp(-e)))

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
      h += j1.T.dot(w[:, None] * j1) / RESOLUTION_DEPENDENT_SIGMA_FACTOR ** 2
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
