import sys

import libtbx.load_env
from libtbx import test_utils

# tst_imports.py only imports modules (no wx.App) and is safe everywhere.
# tst_headless_construct.py and tst_process_control.py create a wx.App, which
# aborts on the no-DISPLAY Linux CI VMs; classify those as expected_unstable
# there so a missing display does not block the tree (mirrors crys3d/run_tests.py).
tst_list = []
tst_list_expected_unstable = []
try:
  import wx  # noqa: F401
  tst_list.append("$D/regression/tst_imports.py")
  _gui_tst_list = [
    "$D/regression/tst_headless_construct.py",
    "$D/regression/tst_process_control.py",
  ]
  if sys.platform in ("darwin", "win32"):
    tst_list.extend(_gui_tst_list)
  else: # no DISPLAY on Azure VMs running linux, so these fail by default
    tst_list_expected_unstable.extend(_gui_tst_list)
except ImportError:
  pass

def run():
  build_dir = libtbx.env.under_build("wxtbx")
  dist_dir = libtbx.env.dist_path("wxtbx")
  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
