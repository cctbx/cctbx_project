import sys

def import_ext(name):
  components = name.split(".")
  if (len(components) > 1):
    __import__(".".join(components[:-1]))
  previous_dlopenflags = None
  if (sys.platform == "linux2"):
    previous_dlopenflags = sys.getdlopenflags()
    sys.setdlopenflags(0x100|0x2)
  mod = __import__(name)
  for comp in components[1:]:
    mod = getattr(mod, comp)
  if (previous_dlopenflags is not None):
    sys.setdlopenflags(previous_dlopenflags)
  return mod

meta_class = import_ext("boost_python_meta_ext").empty.__class__

class injector(object):
  "see boost/libs/python/doc/tutorial/doc/quickstart.txt"

  class __metaclass__(meta_class):

    def __init__(self, name, bases, dict):
      for b in bases:
        if type(b) not in (self, type):
          for k,v in dict.items():
            if (k in ("__init__", "__del__", "__module__", "__file__")):
              continue
            setattr(b,k,v)
      return type.__init__(self, name, bases, dict)
