from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.int_overflow

import boost_adaptbx.boost.python as bp
import sys

def run(args):
  assert len(args) == 0
  for itype in ["int", "long"]:
    sizeof_itype = bp.c_sizeof(itype)
    imax = 2**(8*sizeof_itype-1)-1
    print("Now adding %s values %d + 1 ..." % (itype, imax))
    sys.stdout.flush()
    result = getattr(bp.ext, "add_%ss" % itype)(imax, 1)
    print("Result:", result)

if (__name__ == "__main__"):
  run(sys.argv[1:])
