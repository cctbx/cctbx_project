from __future__ import absolute_import, division, print_function

import boost_adaptbx.python
import sys

def run(args):
  forever = "--forever" in args
  while True:
    boost_adaptbx.python.ext.libtbx_introspection_show_stack()
    boost_adaptbx.python.ext.boost_adaptbx_libc_backtrace(0)
    if (not forever): break
  print("OK")

if (__name__ == "__main__"):
  run(sys.argv[1:])
