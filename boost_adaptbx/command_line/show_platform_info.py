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
  c = getattr(boost.python.ext, "str_or_unicode_as_char_list", None)
  if (c is not None):
    print '"hello" = ', c("hello")
    print 'u"hello" = ', c(u"hello")
    print 'U"hello" = ', c(U"hello")

if (__name__ == "__main__"):
  run()
