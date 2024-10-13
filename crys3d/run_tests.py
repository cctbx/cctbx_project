from __future__ import absolute_import, division, print_function
import sys
from libtbx import test_utils
import libtbx.load_env


tst_list = [ "$D/regression/tst_hklinfo.py" ]
tst_list_expected_unstable = [
   # fails sometimes due to websocket connection problem to webbrowser
   "$D/regression/tst_websocket.py",
   "$D/regression/tst_HKLviewerOSbrowserSliceK-9.py",
   "$D/regression/tst_HKLviewerOSbrowserBinFSigF.py",
]
other_tests = [
  "$D/regression/tst_HKLviewerQtGuiSliceK-9.py",
  "$D/regression/tst_HKLviewerQtGuiBinFSigF.py",
]
try:
  import PySide2  # special import
  pyside2_available = True
except ImportError:
  pyside2_available = False
if (sys.platform == "darwin" or sys.platform == "win32") and pyside2_available:
  tst_list.extend(other_tests)
else: # no DISPLAY environment on Azure VMs running linux so tests fail by default
  tst_list_expected_unstable.extend(other_tests)

# expected failure for Python 2
if sys.version_info < (3, 0):
  tst_list_expected_failures = tst_list
  tst_list_expected_unstable = []
  tst_list = []


def run():
  build_dir = libtbx.env.under_build("crys3d")
  dist_dir = libtbx.env.dist_path("crys3d")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
