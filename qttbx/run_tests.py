from __future__ import absolute_import, division, print_function
from libtbx import test_utils
import libtbx.load_env

tst_list_base = [
  "$D/regression/tst_selection_viewer.py",
  ]
tst_list = tst_list_base

def run():
  build_dir = libtbx.env.under_build("qttbx")
  dist_dir = libtbx.env.dist_path("qttbx")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
