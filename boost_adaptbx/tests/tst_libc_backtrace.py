from __future__ import absolute_import, division, print_function

import boost_adaptbx.boost.python as bp
import sys

def run(args):
  forever = "--forever" in args
  while True:
    bp.ext.libtbx_introspection_show_stack()
    bp.ext.boost_adaptbx_libc_backtrace(0)
    if (not forever): break
  print("OK")

if (__name__ == "__main__"):
  run(sys.argv[1:])
