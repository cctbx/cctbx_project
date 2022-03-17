from __future__ import absolute_import, division, print_function
from libtbx import test_utils
import libtbx.load_env

#tst_list = [
#  "$D/regression/tst_py_from_html.py"
#  ]

tst_list = [
  "$D/regression/tst_1_template.py",
  "$D/regression/tst_2_doc_hlo_intro.py",
  "$D/regression/tst_3_doc_hlo_model_manager.py",
  "$D/regression/tst_4_doc_hlo_data_manager.py",
  "$D/regression/tst_5_doc_hlo_map_manager.py",
  "$D/regression/tst_6_doc_hlo_model_map_manager.py",
  "$D/regression/tst_7_doc_low_flex_advanced.py",
  "$D/regression/tst_8_doc_maps_intro.py",
  "$D/regression/tst_9_doc_maps_boxing.py",
  "$D/regression/tst_10_doc_programming_tips_1.py",
  "$D/regression/tst_11_script_1.py",
  "$D/regression/tst_12_script_compare_ss.py",
  "$D/regression/tst_13_script_ideal_ss.py",
  "$D/regression/tst_14_script_lbfgs_no_curvature.py",
  "$D/regression/tst_15_doc_models_hierarchy.py",
  "$D/regression/tst_16_script_lbfgs_with_curvature.py",
  "$D/regression/tst_17_script_rfactors.py",
  "$D/regression/tst_18_doc_programming_tips_2.py",
  "$D/regression/tst_19_doc_programming_tips_3.py",
  ]

def run():

  build_dir = libtbx.env.under_build("cctbx_website")
  dist_dir = libtbx.env.dist_path("cctbx_website")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
