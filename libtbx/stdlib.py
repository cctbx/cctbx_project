import imp
import sys

def import_module(name):
  try:
    return sys.modules[name]
  except KeyError:
    pass
  if (imp.is_builtin(name)):
    try: return imp.init_builtin(name)
    except: pass
  sys_path = sys.path[:]
  sys_path.reverse()
  for path in sys_path:
    try:
      fp, pathname, description = imp.find_module(name, [path])
    except ImportError:
      pass
    else:
      try:
        return imp.load_module(name, fp, pathname, description)
      finally:
        if (fp is not None): fp.close()
  raise RuntimeError("Cannot import %s module." % name)

math = import_module("math")
random = import_module("random")
