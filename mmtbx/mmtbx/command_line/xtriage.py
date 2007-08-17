# LIBTBX_SET_DISPATCHER_NAME phenix.xtriage

from mmtbx.scaling import xtriage
import sys

if (__name__ == "__main__"):
  xtriage.run(args=sys.argv[1:])
