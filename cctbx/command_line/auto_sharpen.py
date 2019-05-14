from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.auto_sharpen

import sys

if (__name__ == "__main__"):
  from cctbx.maptbx.auto_sharpen import run
  run(args=sys.argv[1:])
