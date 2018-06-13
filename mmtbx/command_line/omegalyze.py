from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.omegalyze
# LIBTBX_SET_DISPATCHER_NAME molprobity.omegalyze

import sys
from mmtbx.validation import omegalyze

if __name__ == "__main__":
  omegalyze.run(sys.argv[1:])
