from __future__ import absolute_import, division, print_function

from libtbx import test_utils
import libtbx.load_env

tst_list = [
  "$B/regression/cpp/tst_infrastructure",
  "$B/regression/cpp/tst_string_view",
  "$B/regression/cpp/tst_tokenizer",
  "$B/regression/cpp/tst_mapped_file",
  "$D/regression/tst_example_tokenize.py",
  "$B/regression/cpp/tst_data_model",
  "$B/regression/cpp/tst_error_handling",
  "$B/regression/cpp/tst_numeric",
  "$D/regression/tst_example_parse.py",
]

def run():
  build_dir = libtbx.env.under_build("xcif")
  dist_dir = libtbx.env.dist_path("xcif")
  test_utils.run_tests(build_dir, dist_dir, tst_list)

if __name__ == '__main__':
  run()
