from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.cablam_training

import sys
from mmtbx.cablam import cablam_training

if __name__ == "__main__":
  cablam_training.run(sys.argv[1:])
