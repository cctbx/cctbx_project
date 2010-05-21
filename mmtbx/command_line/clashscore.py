# LIBTBX_SET_DISPATCHER_NAME phenix.clashscore

import sys
from mmtbx.validation.clashscore import clashscore

if __name__ == "__main__":
  clashscore().run(sys.argv[1:])
