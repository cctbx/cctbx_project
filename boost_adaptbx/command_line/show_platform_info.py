# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.show_platform_info

import boost.python
from libtbx import introspection
import platform
import sys, os

def run():
  sys.stdout.write(boost.python.platform_info)
  print "number of processors:", introspection.number_of_processors(
    return_value_if_unknown="unknown")
  print "os.name:", os.name
  print "sys.platform:", sys.platform
  print "sys.byteorder:", sys.byteorder
  print "platform.platform():", platform.platform()
  print "platform.architecture():", platform.architecture()
  for attr in ["division_by_zero", "invalid", "overflow"]:
    attr = "floating_point_exceptions.%s_trapped" % attr
    print "%s:" % attr, eval("boost.python.%s" % attr)
  try: import thread
  except ImportError: print "import thread: NO"
  else: print "import thread: OK"
  c = getattr(boost.python.ext, "str_or_unicode_as_char_list", None)
  if (c is not None):
    print '"hello" =', c("hello")
    print 'u"hello" =', c(u"hello")
    e = u"\u00C5".encode("utf-8", "strict")
    print 'u"\u00C5" =', c(u"\u00C5"), 'as utf-8 =', c(e)
    print "LATIN CAPITAL LETTER A WITH RING ABOVE =", e

if (__name__ == "__main__"):
  run()
