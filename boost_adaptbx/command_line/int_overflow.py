# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.int_overflow

import boost.python
import sys

def run(args):
  assert len(args) == 0
  sizeof_int = boost.python.c_sizeof("int")
  int_max = 2**(8*sizeof_int-1)-1
  print "Now adding int values %d + 1 ..." % int_max
  sys.stdout.flush()
  result = boost.python.ext.add_ints(int_max, 1)
  print "Result:", result

if (__name__ == "__main__"):
  run(sys.argv[1:])
