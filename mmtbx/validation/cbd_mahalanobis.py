"""
Conformation-dependent Mahalanobis scoring of C-beta deviations.

Each deviation is treated as a 2D vector (x = dev*cos(dihedralNABB),
y = dev*sin(dihedralNABB)) and compared with the distribution of deviation
vectors observed in the same (phi, psi) bin of a high-resolution reference
set. The reference covariances are tighter than the coordinate error of
most models, so an isotropic model-error term is added,

  cov_eff = cov_bin + s**2 * I,

and s is fitted per model so that the median d**2 equals the median of a
chi-squared distribution with two degrees of freedom (2 ln 2).
"""
from __future__ import absolute_import, division, print_function
from math import cos, sin, sqrt, log, radians

from mmtbx.conformation_dependent_library.cdl_utils import round_to_ten
from mmtbx.validation.mean_devs_phi_psi_covariance import cov_devs

median_chi2_2 = 2*log(2)
min_residues_for_calibration = 20

def get_table_class(resname):
  if resname == 'PRO': return 'PRO'
  if resname in ['VAL', 'THR', 'ILE']: return 'VAL_THR_ILE'
  return 'others'

def get_bin_residual(resname, phi_psi, deviation, dihedral):
  """Return (dx, dy, cov_xx, cov_xy, cov_yy) for a residue or None if the
  residue has no phi/psi or no table entry."""
  if not phi_psi or phi_psi[0] is None or phi_psi[1] is None: return None
  key = (round_to_ten(phi_psi[0]), round_to_ten(phi_psi[1]))
  entry = cov_devs[get_table_class(resname)].get(key, None)
  if entry is None: return None
  n, pooled, mu_x, mu_y, c_xx, c_xy, c_yy = entry
  t = radians(dihedral)
  return (deviation*cos(t)-mu_x, deviation*sin(t)-mu_y, c_xx, c_xy, c_yy)

def mahalanobis_d2(residual, s2=0.):
  dx, dy, c_xx, c_xy, c_yy = residual
  c_xx += s2
  c_yy += s2
  det = c_xx*c_yy - c_xy*c_xy
  return (c_yy*dx*dx - 2*c_xy*dx*dy + c_xx*dy*dy)/det

def _median(values):
  v = sorted(values)
  n = len(v)
  if n%2: return v[n//2]
  return (v[n//2-1]+v[n//2])/2

def fit_model_error(residuals, start_s=0.5, cycles=60):
  """Return s (Angstrom) such that the median d**2 is 2 ln 2, or None when
  there are too few residues to calibrate."""
  if len(residuals) < min_residues_for_calibration: return None
  def median_d2(s):
    return _median([mahalanobis_d2(r, s*s) for r in residuals])
  if median_d2(0.) <= median_chi2_2: return 0.
  lo, hi = 0., start_s
  while median_d2(hi) > median_chi2_2: hi *= 2
  for i in range(cycles):
    mid = (lo+hi)/2
    if median_d2(mid) > median_chi2_2: lo = mid
    else: hi = mid
  return (lo+hi)/2

def score(residuals):
  """Return (s, list of d) for a list of residuals (entries may be None).
  When there are too few residues to calibrate, s is None and the
  uncalibrated distances are returned."""
  scored = [r for r in residuals if r is not None]
  s = fit_model_error(scored)
  s2 = (s or 0.)**2
  d = [None if r is None else sqrt(mahalanobis_d2(r, s2)) for r in residuals]
  return s, d
