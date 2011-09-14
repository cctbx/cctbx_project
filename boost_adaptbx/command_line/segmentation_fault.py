# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.segmentation_fault

import boost.python
import sys

def run(args):
  assert len(args) == 0
  print "Now dereferencing null-pointer ..."
  sys.stdout.flush()
  result = boost.python.ext.dereference_char_pointer(None)
  print "Result:", result

if (__name__ == "__main__"):
  run(sys.argv[1:])
