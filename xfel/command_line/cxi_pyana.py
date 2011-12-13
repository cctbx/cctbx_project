# LIBTBX_SET_DISPATCHER_NAME cxi.pyana
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

import os
import sys
from libtbx import easy_run

def run(args):
  sit_reldir = os.environ.get("SIT_RELDIR")
  sit_release = os.environ.get("SIT_RELEASE")
  sit_arch = os.environ.get("SIT_ARCH")
  assert [sit_reldir, sit_release, sit_arch].count(None) == 0
  pyana_path = os.path.join(
    sit_reldir, sit_release, "arch", sit_arch, "bin", "pyana")
  cmd = " ".join(["libtbx.python %s" %pyana_path] + args)
  print cmd
  easy_run.call(cmd)


if __name__ == '__main__':
  run(sys.argv[1:])
