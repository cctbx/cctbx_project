import libtbx.load_env
from libtbx import test_utils

tst_list_base = []
try:
  import PySide2  # noqa: F401
  tst_list_base.append("$D/regression/tst_phil_widgets.py")
except ImportError:
  pass

tst_list = tst_list_base

def run():
  build_dir = libtbx.env.under_build("qttbx")
  dist_dir = libtbx.env.dist_path("qttbx")
  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
