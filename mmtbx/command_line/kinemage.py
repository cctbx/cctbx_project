# LIBTBX_SET_DISPATCHER_NAME phenix.kinemage

import sys
from mmtbx.kinemage import validation

if __name__ == "__main__":
  validation.run(sys.argv[1:])
