from libtbx import test_utils
import libtbx.load_env

def run():
  tst_list = (
  "$D/test_utils.py",
  "$D/utils.py",
  "$D/sge_utils.py",
  "$D/introspection.py",
  "$D/easy_run.py",
  "$D/tst_utils.py",
  "$D/tst_math_utils.py",
  "$D/assert_utils.py",
  "$D/str_utils.py",
  "$D/table_utils.py",
  "$D/tst_dlite.py",
  "$D/phil/tst_tokenizer.py",
  "$D/phil/tst.py",
  "$D/tst_object_oriented_patterns.py",
  "$D/tst_itertbx.py",
  "$D/tst_symmetric_multi_processing.py",
  )

  build_dir = None
  dist_dir = libtbx.env.dist_path("libtbx")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
