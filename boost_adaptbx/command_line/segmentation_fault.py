from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.segmentation_fault

import boost_adaptbx.boost.python as bp
import sys

def run(args):
  assert len(args) == 0
  print("Now dereferencing null-pointer ...")
  sys.stdout.flush()
  result = bp.ext.dereference_char_pointer(None)
  print("Result:", result)

if (__name__ == "__main__"):
  run(sys.argv[1:])
