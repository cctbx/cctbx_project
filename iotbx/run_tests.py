import os, os.path
from scitbx import test_utils

def run():
  tst_list = (
  "$D/iotbx/pdb/tst_pdb.py",
  "$D/iotbx/mtz/tst_mtz.py",
  ["$D/iotbx/mtz/regression.py", "P31"],
  "$D/iotbx/detectors/tst_adsc.py",
  "$D/iotbx/xplor/tst_xplormap.py",
  )

  build_dir = os.path.join(os.environ["LIBTBX_BUILD"], "iotbx")
  dist_dir = os.environ["IOTBX_DIST"]

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
