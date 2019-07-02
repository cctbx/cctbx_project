from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.fobs_minus_fobs_map

import mmtbx.maps.fobs_minus_fobs_map
import sys

if(__name__ == "__main__"):
  mmtbx.maps.fobs_minus_fobs_map.run(args = sys.argv[1:])
