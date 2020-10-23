from __future__ import absolute_import, division, print_function
from libtbx import test_utils
import libtbx.load_env

#tst_list = [
#  "$D/regression/tst_py_from_html.py"
#  ]

tst_list = [
  "$D/regression/tst_1_template.py",
  "$D/regression/tst_2_doc_high_level_objects.py",
  "$D/regression/tst_3_doc_model_manager.py",
  "$D/regression/tst_4_doc_data_manager.py",
  "$D/regression/tst_5_doc_map_manager.py",
  "$D/regression/tst_6_doc_model_map_manager.py",
  ]

def run():

  build_dir = libtbx.env.under_build("cctbx_website")
  dist_dir = libtbx.env.dist_path("cctbx_website")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
