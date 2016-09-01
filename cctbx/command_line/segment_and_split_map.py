from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.segment_and_split_map

import sys

if (__name__ == "__main__"):
  from cctbx.maptbx.segment_and_split_map import run
  run(args=sys.argv[1:])
