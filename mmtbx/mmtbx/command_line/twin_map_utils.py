# LIBTBX_SET_DISPATCHER_NAME phenix.twin_map_utils

from mmtbx.twinning import twin_map_utils
import sys

if (__name__ == "__main__"):
  twin_map_utils.run(args=sys.argv[1:])
