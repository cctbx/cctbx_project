from __future__ import generators
import inspect

__builtins__.setdefault("False", 0)
__builtins__.setdefault("True", 1)

if (__builtins__.get("bool", None) is None):
  def bool(value):
    if (value): return True
    return False
  __builtins__["bool"] = bool

if (__builtins__.get("_enumerate", None) is None):
  def enumerate(iterable):
    """enumerate(iterable) -> iterator for index, value of iterable

  Return an enumerate object.  iterable must be an other object that supports
  iteration.  The enumerate object yields pairs containing a count (from
  zero) and a value yielded by the iterable argument.  enumerate is useful
  for obtaining an indexed list: (0, seq[0]), (1, seq[1]), (2, seq[2]), ..."""
    i = 0
    it = iter(iterable)
    while 1:
      yield (i, it.next())
      i += 1
  __builtins__["enumerate"] = enumerate

# libtbx enhancement
if (__builtins__.get("signatures", None) is None):
  def signatures(function_object):
    signatures = getattr(function_object, "__signatures__", None)
    if (signatures is None):
      im_func = getattr(function_object, "im_func", None)
      if (im_func is not None):
        signatures = getattr(im_func, "__signatures__", None)
    if (signatures is not None):
      return signatures
    name = []
    try: module = inspect.getmodule(function_object)
    except KeyboardInterrupt: raise
    except: pass
    else:
      if (module is not None):
        name.append(module.__name__)
    name.append(function_object.__name__)
    return [".".join(name) + inspect.formatargspec(
      *inspect.getargspec(function_object))]
  __builtins__["signatures"] = signatures
