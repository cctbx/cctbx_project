import os, os.path
from libtbx import test_utils

def run():
  tst_list = (
  "$D/iotbx/tst_wildcard.py",
  "$D/iotbx/tst_simple_parser.py",
  "$D/iotbx/tst_phil.py",
  "$D/iotbx/kriber/tst_strudat.py",
  "$D/iotbx/pdb/tst_pdb.py",
  "$D/iotbx/cns/space_group_symbols.py",
  "$D/iotbx/cns/tst_cns.py",
  ["$D/iotbx/scalepack/tst_merge.py", "P31"],
  "$D/include/iotbx/mtz/tst_ext.py",
  "$D/iotbx/mtz/extract_from_symop_lib.py",
  ["$D/iotbx/mtz/tst.py", "P31"],
  "$D/iotbx/tst_reflection_file_utils.py",
  "$D/iotbx/detectors/tst_adsc.py",
  "$D/iotbx/xplor/tst_xplormap.py",
  ["$D/iotbx/regression/tst_reflection_statistics.py", "P31m"]
  )

  build_dir = os.path.join(os.environ["LIBTBX_BUILD"], "iotbx")
  dist_dir = os.environ["IOTBX_DIST"]

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
