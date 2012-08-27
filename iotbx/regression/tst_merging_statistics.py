
from __future__ import division
from libtbx import easy_run
import libtbx.load_env
import os
import sys

def exercise (debug=False) :
  if (not libtbx.env.has_module("phenix_regression")) :
    print "phenix_regression not configured, skipping."
    return
  hkl_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/wizards/p9_se_w2.sca",
    test=os.path.isfile)
  args = [
    "phenix.merging_statistics",
    hkl_file,
    "space_group=I4",
    "unit_cell=113.949,113.949,32.474,90,90,90",
  ]
  if (debug) :
    args.append("debug=True")
    print " ".join(args)
  result = easy_run.fully_buffered(" ".join(args)).raise_if_errors()
  assert ("R-merge: 0.073" in result.stdout_lines)
  assert ("R-meas:  0.079" in result.stdout_lines)
  print "OK"

if (__name__ == "__main__") :
  exercise(debug=("--debug" in sys.argv))
