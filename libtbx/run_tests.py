from libtbx import test_utils
import os, os.path

def run():
  tst_list = (
  "$D/libtbx/str_utils.py",
  "$D/phil/tst.py",
  )

  build_dir = None
  dist_dir = os.environ["LIBTBX_DIST"]

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
