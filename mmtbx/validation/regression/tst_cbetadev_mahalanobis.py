from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
import iotbx.pdb
from mmtbx.validation import cbetadev
from mmtbx.validation import cbd_mahalanobis

def get_hierarchy():
  import mmtbx
  from pathlib import Path
  data_dir = Path(mmtbx.__file__).parent / 'regression' / 'pdbs'
  return iotbx.pdb.input(str(data_dir / '2ERL_noH.pdb')).construct_hierarchy()

def tst_01():
  """
  Scoring functions.
  """
  # no phi/psi or no table entry
  assert cbd_mahalanobis.get_bin_residual('ALA', None, 0.1, 0.) is None
  assert cbd_mahalanobis.get_bin_residual('ALA', (None, 120.), 0.1, 0.) is None
  # deviation equal to the bin mean has d = 0
  key = (-60, -40)
  n, pooled, mu_x, mu_y, c_xx, c_xy, c_yy = cbd_mahalanobis.cov_devs['others'][key]
  import math
  dev = math.sqrt(mu_x**2+mu_y**2)
  dihedral = math.degrees(math.atan2(mu_y, mu_x))
  residual = cbd_mahalanobis.get_bin_residual('LEU', key, dev, dihedral)
  assert approx_equal(cbd_mahalanobis.mahalanobis_d2(residual), 0)
  # identity covariance
  assert approx_equal(cbd_mahalanobis.mahalanobis_d2((3., 4., 1., 0., 1.)), 25)
  assert approx_equal(cbd_mahalanobis.mahalanobis_d2((3., 4., 1., 0., 1.), 1.), 12.5)
  # too few residues to calibrate
  s, ds = cbd_mahalanobis.score([(1., 0., 1., 0., 1.)]*5 + [None])
  assert s is None
  assert ds[-1] is None
  assert approx_equal(ds[0], 1)
  # median d**2 is 2 ln 2 after calibration
  residuals = [(0.1*i, 0., 1e-4, 0., 1e-4) for i in range(1, 30)]
  s, ds = cbd_mahalanobis.score(residuals)
  assert s > 0
  d2 = sorted(d*d for d in ds)
  assert approx_equal(d2[len(d2)//2], 2*math.log(2), eps=1e-6)

def tst_02():
  """
  Mahalanobis scoring is reported alongside the 0.25 Angstrom outliers.
  """
  hierarchy = get_hierarchy()
  default = cbetadev.cbetadev(pdb_hierarchy=hierarchy, out=null_out())
  assert default.mahalanobis is None
  assert default.get_mahalanobis_outlier_count() is None
  for result in default.results:
    assert result.mahalanobis_d is None
  validation = cbetadev.cbetadev(pdb_hierarchy=hierarchy,
                                 mahalanobis_z_score=True,
                                 mahalanobis_cutoff=2.2,
                                 out=null_out())
  assert validation.n_outliers == default.n_outliers == 3
  assert len(validation.results) == len(default.results)
  m = validation.mahalanobis
  assert approx_equal(m.model_error_by_model[''], 0.0726, eps=1e-3)
  assert m.n_scored_by_model[''] == 40
  assert validation.get_mahalanobis_outlier_count() == 4
  outliers = [r for r in validation.results
              if r.mahalanobis_d is not None and r.mahalanobis_d >= 2.2]
  assert len(outliers) == 4
  # Tyr 39 is only a Mahalanobis outlier
  assert len([r for r in outliers if not r.is_outlier()]) == 1
  # outliers_only keeps both kinds of outlier
  outliers_only = cbetadev.cbetadev(pdb_hierarchy=hierarchy,
                                    outliers_only=True,
                                    mahalanobis_z_score=True,
                                    mahalanobis_cutoff=2.2,
                                    out=null_out())
  ids = set(r.id_str() for r in outliers_only.results)
  expected = set(r.id_str() for r in validation.results
                 if r.is_outlier() or r in outliers)
  assert ids == expected, (ids, expected)

if (__name__ == "__main__"):
  tst_01()
  tst_02()
  print("OK")
