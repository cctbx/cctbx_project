from libtbx import test_utils
import os

def run():
  tst_list = (
  "$D/tst_rational.py",
  )

  build_dir = os.path.join(os.environ["LIBTBX_BUILD"], "boost_adaptbx")
  dist_dir = os.environ["BOOST_ADAPTBX_DIST"]

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
