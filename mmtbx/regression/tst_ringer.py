
from __future__ import division
from libtbx import easy_run
import libtbx.load_env # import dependency
import os

def exercise () :
  if (not libtbx.env.has_module("phenix_regression")) :
    print "phenix_regression not available, skipping test."
    return
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1yjp_h.pdb",
    test=os.path.isfile)
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/1yjp.mtz",
    test=os.path.isfile)
  result = easy_run.fully_buffered(
    "phenix.maps \"%s\" \"%s\" output.prefix=1yjp" %
    (pdb_file, mtz_file)).raise_if_errors()
  assert (result.return_code == 0)
  result = easy_run.fully_buffered(
    "mmtbx.ringer \"%s\" 1yjp_map_coeffs.mtz" % pdb_file).raise_if_errors()
  assert (result.return_code == 0)
  print "OK"

if (__name__ == "__main__") :
  exercise()
