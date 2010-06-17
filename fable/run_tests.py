from libtbx import test_utils
import libtbx.load_env

def run():
  tst_list = (
  "$D/tst_ext.py",
  "$D/tst_equivalence.py",
  "$D/tst_read.py",
  "$D/tst_cout.py",
  "$D/test/tst_show_calls.py",
  "$D/test/tst_command_line.py",
  "$D/test/tst_separate_files.py",
  "$D/test/tst_io.py",
  ["$D/tst_cout_compile.py", "stop"],
  )

  build_dir = libtbx.env.under_build("fable")
  dist_dir = libtbx.env.dist_path("fable")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
