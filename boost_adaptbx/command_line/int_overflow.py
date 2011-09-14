# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.int_overflow

import boost.python
import sys

def run(args):
  assert len(args) == 0
  for itype in ["int", "long"]:
    sizeof_itype = boost.python.c_sizeof(itype)
    imax = 2**(8*sizeof_itype-1)-1
    print "Now adding %s values %d + 1 ..." % (itype, imax)
    sys.stdout.flush()
    result = getattr(boost.python.ext, "add_%ss" % itype)(imax, 1)
    print "Result:", result

if (__name__ == "__main__"):
  run(sys.argv[1:])
