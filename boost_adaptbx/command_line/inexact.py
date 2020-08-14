from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.inexact

import boost_adaptbx.boost.python as bp
import sys

def run(args):
  assert len(args) == 0
  print("Now creating a NaN in C++ as 0/0 ...")
  sys.stdout.flush()
  result = bp.ext.divide_doubles(0, 0)
  print("Result:", result)

if (__name__ == "__main__"):
  run(sys.argv[1:])
