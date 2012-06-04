from libtbx import test_utils
import libtbx.load_env

tst_list = (
  "$B/tests/tst_optional_copy",
  "$D/tests/tst_libc_backtrace.py",
  "$D/tests/tst_rational.py",
  "$D/tests/tst_rational_truediv.py",
  "$D/tests/tst_optional.py",
  "$D/tests/tst_std_pair.py",
  "$D/tests/tst_tuple.py",
  "$D/tests/tst_python_streambuf.py",
  "$D/tests/tst_stdout.py",
  "$D/tests/tst_stderr_stdout.py",
  )


def run():
  build_dir = libtbx.env.under_build("boost_adaptbx")
  dist_dir = libtbx.env.dist_path("boost_adaptbx")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
