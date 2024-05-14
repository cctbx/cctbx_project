
from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import libtbx.load_env
import os

def exercise_simple():
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/ha_patterson.mtz",
    test=os.path.isfile)
  if (mtz_file is None):
    print("phenix_regression not available, skipping")
    return
  result = easy_run.fully_buffered("cctbx.patterson_map \"%s\"" % mtz_file
    ).raise_if_errors()
  assert (result.return_code == 0)

if (__name__ == "__main__"):
  exercise_simple()
  print("OK")
