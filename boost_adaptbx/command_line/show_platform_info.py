# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.show_platform_info

import boost.python
from libtbx import introspection
import platform
import sys, os

def run (out=None, omit_unicode_experiment=False):
  if (out is None): out = sys.stdout
  out.write(boost.python.platform_info)
  print >> out, "os.name:", os.name
  print >> out, "sys.platform:", sys.platform
  print >> out, "sys.byteorder:", sys.byteorder
  print >> out, "platform.platform():", platform.platform()
  print >> out, "platform.architecture():", platform.architecture()
  for attr in ["division_by_zero", "invalid", "overflow"]:
    attr = "floating_point_exceptions.%s_trapped" % attr
    print >> out, "%s:" % attr, eval("boost.python.%s" % attr)
  print >> out, "number of processors:", introspection.number_of_processors(
    return_value_if_unknown="unknown")
  introspection.machine_memory_info().show(out=out)
  try: import thread
  except ImportError: print >> out, "import thread: NO"
  else: print >> out, "import thread: OK"
  print "Division operator semantics: %s division" \
    % ["floor", "true"][int(1/2 != 0)]
  c = getattr(boost.python.ext, "str_or_unicode_as_char_list", None)
  if (c is not None and not omit_unicode_experiment):
    print >> out, '"hello" =', c("hello")
    print >> out, 'u"hello" =', c(u"hello")
    e = u"\u00C5".encode("utf-8", "strict")
    print >> out, 'u"\u00C5" =', c(u"\u00C5"), 'as utf-8 =', c(e)
    print >> out, "LATIN CAPITAL LETTER A WITH RING ABOVE =", e

if (__name__ == "__main__"):
  run()
