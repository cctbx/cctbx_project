from libtbx import test_utils
import libtbx.load_env

tst_list = (
  "$D/core_toolbox/boost_python/tst_small_image.py",
  "$D/diffraction/geometry.py",
  "$D/math_support/sphere_formulae.py"
  )

def run_standalones():
  build_dir = libtbx.env.under_build("spotfinder")
  dist_dir = libtbx.env.dist_path("spotfinder")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run_standalones()
