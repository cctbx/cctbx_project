import os, os.path
from scitbx import test_utils

def run():
  tst_list = (
  "$B/array_family/tst_af_1.exe",
  "$B/array_family/tst_af_2.exe",
  "$B/array_family/tst_af_3.exe",
  "$B/array_family/tst_af_4.exe",
  "$B/array_family/tst_af_5.exe",
  "$B/array_family/tst_vec3.exe",
  "$B/array_family/tst_mat3.exe",
  "$B/array_family/tst_sym_mat3.exe",
  "$B/array_family/tst_mat_ref.exe",
  "$B/array_family/tst_accessors.exe",
  "$D/array_family/boost_python/regression_test.py",
  "$D/array_family/boost_python/tst_flex.py",
  "$D/boost_python/tst_rational.py",
  "$B/serialization/tst_base_256.exe",
  ["$D/lbfgs/boost_python/tst_lbfgs.py"],
  ["$D/fftpack/boost_python/tst_fftpack.py"],
  )

  build_dir = os.path.join(os.environ["LIBTBX_BUILD"], "scitbx")
  dist_dir = os.environ["SCITBX_DIST"]

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
