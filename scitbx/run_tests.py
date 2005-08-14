from libtbx import test_utils
import libtbx.load_env

def run():
  tst_list = (
  "$B/array_family/tst_af_1",
  "$B/array_family/tst_af_2",
  "$B/array_family/tst_af_3",
  "$B/array_family/tst_af_4",
  "$B/array_family/tst_af_5",
  "$B/array_family/tst_vec3",
  "$B/array_family/tst_mat3",
  "$B/array_family/tst_sym_mat3",
  "$B/array_family/tst_mat_ref",
  "$B/array_family/tst_accessors",
  "$B/serialization/tst_base_256",
  "$D/include/scitbx/stl/tst_map.py",
  "$D/include/scitbx/stl/tst_set.py",
  "$D/include/scitbx/stl/tst_vector.py",
  "$D/array_family/boost_python/regression_test.py",
  "$D/array_family/boost_python/tst_flex.py",
  "$D/array_family/boost_python/tst_shared.py",
  "$D/scitbx/matrix.py",
  "$D/scitbx/python_utils/math_utils.py",
  "$D/scitbx/python_utils/tst_random_transform.py",
  "$D/math/boost_python/tst_math.py",
  "$D/scitbx/math/tst_superpose.py",
  ["$D/lbfgs/boost_python/tst_lbfgs.py"],
  "$D/lbfgsb/boost_python/tst_lbfgsb.py",
  ["$D/fftpack/boost_python/tst_fftpack.py"],
  "$D/scitbx/examples/lbfgs_recipe.py",
  "$D/scitbx/examples/lbfgs_linear_least_squares_fit.py",
  "$D/scitbx/examples/chebyshev_lsq_example.py"
  )

  build_dir = libtbx.env.under_build("scitbx")
  dist_dir = libtbx.env.dist_path("scitbx")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
