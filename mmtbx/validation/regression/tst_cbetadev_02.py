
from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
import libtbx.load_env
import iotbx.pdb
from mmtbx.validation import cbetadev
import os

def tst_01():
  """
  Exercise with model where coordinates are rounded to whole angstroms.
  """
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="cctbx_project/mmtbx/regression/pdbs/1ucs_cutted_xyz_rounded.pdb",
    test=os.path.isfile)

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
