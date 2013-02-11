
from __future__ import division
from libtbx import easy_run
import libtbx.load_env
import os

def exercise () :
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/ha_patterson.mtz",
    test=os.path.isfile)
  if (mtz_file is None) :
    print "phenix_regression not available, skipping"
    return
  result = easy_run.fully_buffered("cctbx.patterson_map \"%s\"" % mtz_file
    ).raise_if_errors()
  assert (result.return_code == 0)
  print "OK"

if (__name__ == "__main__") :
  exercise()
