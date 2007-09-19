# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.show_platform_info

import boost.python
import sys

def run():
  sys.stdout.write(boost.python.platform_info)
  print "sys.byteorder:", sys.byteorder
  try: import thread
  except ImportError: print "import thread: NO"
  else: print "import thread: OK"

if (__name__ == "__main__"):
  run()
