from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.omegalyze

import sys
from mmtbx.validation import omegalyze

if __name__ == "__main__":
  omegalyze.run(sys.argv[1:])
