import os, os.path
from libtbx import test_utils

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
  "$D/array_family/boost_python/regression_test.py",
  "$D/array_family/boost_python/tst_flex.py",
  "$B/serialization/tst_base_256",
  "$D/scitbx/python_utils/math_utils.py",
  "$D/math/boost_python/tst_math.py",
  ["$D/lbfgs/boost_python/tst_lbfgs.py"],
  "$D/lbfgsb/boost_python/tst_lbfgsb.py",
  ["$D/fftpack/boost_python/tst_fftpack.py"],
  )

  build_dir = os.path.join(os.environ["LIBTBX_BUILD"], "scitbx")
  dist_dir = os.environ["SCITBX_DIST"]

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
