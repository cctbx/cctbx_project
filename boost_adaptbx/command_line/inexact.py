# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.inexact

import boost.python
import sys

def run(args):
  assert len(args) == 0
  print "Now creating a NaN in C++ as 0/0 ..."
  sys.stdout.flush()
  result = boost.python.ext.divide_doubles(0, 0)
  print "Result:", result

if (__name__ == "__main__"):
  run(sys.argv[1:])
