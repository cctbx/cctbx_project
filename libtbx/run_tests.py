from libtbx import test_utils
import libtbx.load_env

def run():
  tst_list = (
  "$D/libtbx/test_utils.py",
  "$D/libtbx/utils.py",
  "$D/libtbx/sge_utils.py",
  "$D/libtbx/introspection.py",
  "$D/libtbx/easy_run.py",
  "$D/libtbx/tst_utils.py",
  "$D/libtbx/tst_math_utils.py",
  "$D/libtbx/assert_utils.py",
  "$D/libtbx/str_utils.py",
  "$D/libtbx/table_utils.py",
  "$D/libtbx/tst_dlite.py",
  "$D/libtbx/phil/tst_tokenizer.py",
  "$D/libtbx/phil/tst.py",
  "$D/libtbx/tst_object_oriented_patterns.py",
  "$D/libtbx/tst_itertbx.py",
  "$D/libtbx/tst_symmetric_multi_processing.py",
  )

  build_dir = None
  dist_dir = libtbx.env.dist_path("libtbx")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
