from __future__ import division
from libtbx import test_utils
import libtbx.load_env

tst_list = ["$D/test/command_line/tst_cxi_index.py"]

def run():
  build_dir = libtbx.env.under_build("xfel")
  dist_dir = libtbx.env.dist_path("xfel")
  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()

# libtbx.run_tests_parallel module=xfel
