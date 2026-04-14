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
  "$D/regression/tst_bindings.py",
  "$D/regression/tst_memory.py",
  "$D/regression/tst_gil_release.py",
  "$D/regression/tst_string_interning.py",
  "$D/regression/tst_api_compat.py",
  "$D/regression/tst_parse_verify.py",
  "$D/regression/tst_utf8.py",
  "$D/regression/tst_roundtrip.py",
  "$D/regression/tst_pdb_coordinates.py",
  "$D/regression/tst_pdb_sf.py",
  "$B/regression/cpp/tst_non_strict_mode",
  "$D/regression/tst_non_strict_mode.py",
  "$B/regression/cpp/tst_thread_safety",
  "$B/regression/cpp/fuzz_target",
]

def run():
  build_dir = libtbx.env.under_build("xcif")
  dist_dir = libtbx.env.dist_path("xcif")
  test_utils.run_tests(build_dir, dist_dir, tst_list)

if __name__ == '__main__':
  run()
