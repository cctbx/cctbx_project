from libtbx import test_utils
import libtbx.load_env

def run():
  tst_list = (
  "$D/tst_rational.py",
  "$D/tst_rational_truediv.py",
  "$D/tst_optional.py",
  )

  build_dir = libtbx.env.under_build("boost_adaptbx")
  dist_dir = libtbx.env.dist_path("boost_adaptbx")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
