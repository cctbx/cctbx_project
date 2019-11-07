from __future__ import absolute_import, division, print_function
from libtbx import test_utils
import libtbx.load_env

tst_list = (
  "$D/tst_ext.py",
  "$D/tst_equivalence.py",
  "$D/tst_read.py",
  "$D/tst_cout.py",
  "$D/test/tst_show_calls.py",
  "$D/test/tst_command_line.py",
  "$D/test/tst_separate_files.py",
  ["$D/tst_cout_compile.py", "stop"],
  )

tst_list_expected_unstable = [
  "$D/test/tst_io.py",
]

def run():
  build_dir = libtbx.env.under_build("fable")
  dist_dir = libtbx.env.dist_path("fable")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
