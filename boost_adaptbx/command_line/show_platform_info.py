from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.show_platform_info

import boost_adaptbx.boost.python as bp
from libtbx import introspection
import libtbx.load_env
import platform
import sys, os
op = os.path

class div_probe(object):
  def __div__(O, other): return "floor"
  def __truediv__(O, other): return "true"

def run(out=None, omit_unicode_experiment=False):
  if (out is None): out = sys.stdout
  out.write(bp.platform_info)
  tag = libtbx.env.under_dist("boost", "TAG")
  if (op.isfile(tag)):
    tag = open(tag).read().strip()
  else:
    tag = None
  print("boost/TAG:", tag, file=out)
  print("os.name:", os.name, file=out)
  print("sys.platform:", sys.platform, file=out)
  print("sys.byteorder:", sys.byteorder, file=out)
  print("platform.platform():", platform.platform(), file=out)
  print("platform.architecture():", platform.architecture(), file=out)
  for attr in ["division_by_zero", "invalid", "overflow"]:
    attr = "floating_point_exceptions.%s_trapped" % attr
    print("%s:" % attr, eval("bp.%s" % attr), file=out)
  print("number of processors:", introspection.number_of_processors(
    return_value_if_unknown="unknown"), file=out)
  introspection.machine_memory_info().show(out=out)
  try: import thread
  except ImportError: print("import thread: NO", file=out)
  else: print("import thread: OK", file=out)
  print("Division operator semantics: %s division" % (div_probe() / 0))
  c = getattr(bp.ext, "str_or_unicode_as_char_list", None)
  if (c is not None and not omit_unicode_experiment):
    print('"hello" =', c("hello"), file=out)
    print('u"hello" =', c(u"hello"), file=out)
    e = u"\u00C5".encode("utf-8", "strict")
    # XXX temp fix for python3 failure: unicode failure in show_platform_info.py
    # Reason: c(u"\u00C5") fails and c(e) fails
    try:
      cc = c(u"\u00C5")
      ce = c(e)
    except Exception as ignore_exception:
      cc = u"\u00C5"
      ce = e

    print('u"\u00C5" =', cc, 'as utf-8 =', ce, file=out)
    print("LATIN CAPITAL LETTER A WITH RING ABOVE =", e, file=out)
  from libtbx.utils import format_cpu_times
  print(format_cpu_times())

if (__name__ == "__main__"):
  run()
