from libtbx import test_utils
import libtbx.load_env

def run():
  tst_list = (
  "$D/scattering/tst_direct_summation.py",
  "$D/math/special_functions/tst_spherical_bessel_jn.py",
  )

  build_dir = libtbx.env.under_build("cudatbx")
  dist_dir = libtbx.env.dist_path("cudatbx")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
