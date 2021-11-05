from __future__ import absolute_import, division, print_function
import sys
from libtbx import test_utils
import libtbx.load_env

tst_list = [
"$D/regression/tst_hklinfo.py",
]

# expected failure for Python 2
if sys.version_info < (3, 0):
  tst_list_expected_failures = tst_list
  tst_list = []

def run():
  build_dir = libtbx.env.under_build("crys3d")
  dist_dir = libtbx.env.dist_path("crys3d")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
