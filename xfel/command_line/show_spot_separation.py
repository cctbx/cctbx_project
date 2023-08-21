from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.show_spot_separation

from xfel.util.show_spot_separation import run
import sys

if __name__ == '__main__':
  run(sys.argv[1:])
