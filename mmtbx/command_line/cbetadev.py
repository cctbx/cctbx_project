# LIBTBX_SET_DISPATCHER_NAME phenix.cbetadev

import sys
from mmtbx.validation.cbetadev import cbetadev

if __name__ == "__main__":
  cbetadev().run(sys.argv[1:])
