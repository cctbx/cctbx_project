
from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
import iotbx.pdb
from mmtbx.validation import cbetadev

def tst_01():
  """
  Exercise with model where coordinates are rounded to whole angstroms.
  """
  import mmtbx
  from pathlib import Path
  data_dir = Path(mmtbx.__file__).parent / 'regression' / 'pdbs'
  regression_pdb = str( data_dir / '1ucs_cutted_xyz_rounded.pdb')

  pdb_in = iotbx.pdb.input(regression_pdb)
  hierarchy = pdb_in.construct_hierarchy()
  validation = cbetadev.cbetadev(
    pdb_hierarchy=hierarchy,
    outliers_only=True)
  # print (validation.get_weighted_outlier_percent())
  assert approx_equal(validation.get_weighted_outlier_percent(), 96.7741935483871)

if (__name__ == "__main__"):
  tst_01()
  print("OK")
