from __future__ import absolute_import, division, print_function
import sys
from libtbx import test_utils
import libtbx.load_env

if "darwin" in sys.platform:
  tst_list = [
    "$D/regression/tst_hklinfo.py",
    "$D/regression/tst_websocket.py",
    "$D/regression/tst_HKLviewerOSbrowserSliceK-9.py",
    "$D/regression/tst_HKLviewerOSbrowserBinFSigF.py",
    "$D/regression/tst_HKLviewerQtGuiSliceK-9.py",
    "$D/regression/tst_HKLviewerQtGuiBinFSigF.py",
  ]


if sys.platform == "win32":
  tst_list = [
    "$D/regression/tst_hklinfo.py",
  ]
  tst_list_expected_unstable = [ # browser with websocket has latency issue in a virtual machine
    "$D/regression/tst_websocket.py",
    "$D/regression/tst_HKLviewerOSbrowserSliceK-9.py",
    "$D/regression/tst_HKLviewerOSbrowserBinFSigF.py",
  ]

  tst_list_expected_failures = [ # WebGL + QWebEngine doesn't quite work in a virtual machine
    "$D/regression/tst_HKLviewerQtGuiSliceK-9.py",
    "$D/regression/tst_HKLviewerQtGuiBinFSigF.py",
  ]


if "linux" in sys.platform:
  tst_list = [
    "$D/regression/tst_hklinfo.py",
  ]
  tst_list_expected_unstable = []
  tst_list_expected_failures = [ # no DISPLAY environment on azure
    "$D/regression/tst_websocket.py",
    "$D/regression/tst_HKLviewerOSbrowserSliceK-9.py",
    "$D/regression/tst_HKLviewerOSbrowserBinFSigF.py",
    "$D/regression/tst_HKLviewerQtGuiSliceK-9.py",
    "$D/regression/tst_HKLviewerQtGuiBinFSigF.py",
  ]

# expected failure for Python 2
if sys.version_info < (3, 0):
  tst_list_expected_failures = tst_list
  tst_list = []


def run():
  build_dir = libtbx.env.under_build("crys3d")
  dist_dir = libtbx.env.dist_path("crys3d")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
