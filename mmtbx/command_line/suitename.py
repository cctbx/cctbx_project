# LIBTBX_SET_DISPATCHER_NAME phenix.suitename
# LIBTBX_SET_DISPATCHER_NAME molprobity.suitename
# LIBTBX_SET_DISPATCHER_NAME cctbx.suitename

from __future__ import absolute_import, division, print_function
import sys
from  mmtbx.programs import suitename

if __name__ == '__main__':
  suitename.run(args=sys.argv[1:])

