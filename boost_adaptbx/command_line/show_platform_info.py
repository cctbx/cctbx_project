# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.show_platform_info

import boost.python
from libtbx import introspection
import sys

def run():
  sys.stdout.write(boost.python.platform_info)
  print "number of processors:", introspection.number_of_processors(
    return_value_if_unknown="unknown")
  print "sys.byteorder:", sys.byteorder
  try: import thread
  except ImportError: print "import thread: NO"
  else: print "import thread: OK"

if (__name__ == "__main__"):
  run()
