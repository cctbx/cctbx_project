from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.cablam

import sys
from mmtbx.validation import cablam

if __name__ == "__main__":
  cablam.run(sys.argv[1:])
