# LIBTBX_SET_DISPATCHER_NAME phenix.suitename_old
# LIBTBX_SET_DISPATCHER_NAME molprobity.suitename_old
# LIBTBX_SET_DISPATCHER_NAME cctbx.suitename_old

from __future__ import absolute_import, division, print_function
import sys
from  mmtbx.programs import suitename_old

if __name__ == '__main__':
  suitename_old.run(args=sys.argv[1:])

