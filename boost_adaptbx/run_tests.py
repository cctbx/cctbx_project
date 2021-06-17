from __future__ import absolute_import, division, print_function
from libtbx import test_utils
from libtbx.env_config import get_gcc_version
import sys
import libtbx.load_env

tst_list_base = [
  "$B/tests/tst_optional_copy",
  "$D/tests/tst_libc_backtrace.py",
  "$D/tests/tst_rational.py",
  "$D/tests/tst_optional.py",
  "$D/tests/tst_std_pair.py",
  "$D/tests/tst_tuple.py",
  "$D/tests/tst_stdout.py",
  "$D/tests/tst_stderr_stdout.py",
  "$D/tests/tst_deprecation_tools.py",
  ]

# failing test on Python 3 or GCC > 7
tst_list_fail_py3 = [
  "$D/tests/tst_python_streambuf.py",
  ]
tst_list_fail = list()
gcc_version = get_gcc_version()
if gcc_version is None:
  gcc_version = 0
if sys.version_info[0] > 2 \
   or (sys.platform.startswith('linux') and gcc_version >= 70000):
  tst_list_fail += tst_list_fail_py3
else:
  tst_list_base += tst_list_fail_py3

# final lists
tst_list = tst_list_base
tst_list_expected_failures = tst_list_fail

def run():
  build_dir = libtbx.env.under_build("boost_adaptbx")
  dist_dir = libtbx.env.dist_path("boost_adaptbx")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
