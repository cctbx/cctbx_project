import boost.python
import sys

def run(args):
  forever = "--forever" in args
  while True:
    boost.python.ext.libtbx_introspection_show_stack()
    boost.python.ext.boost_adaptbx_libc_backtrace(0)
    if (not forever): break
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
