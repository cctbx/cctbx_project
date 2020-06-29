from __future__ import absolute_import, division, print_function
from libtbx import test_utils
import libtbx.load_env

tst_list = [
  "$B/tests/tst_fast_linalg",
  ]

def run():
  build_dir = libtbx.env.under_build("fast_linalg")
  dist_dir = libtbx.env.dist_path("fast_linalg")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
