import sys

def adopt_init_args(obj, args, exclude=()):
  del args["self"]
  for param in exclude:
    del args[param]
  for key in args.keys():
    assert not hasattr(obj.__dict__, key)
  obj.__dict__.update(args)

def import_regular_symbols(dict_target, dict_source):
  for key, value in dict_source.items():
    if (not key.startswith("_") and not key in dict_target):
      dict_target[key] = value

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
  if (previous_dlopenflags != None):
    sys.setdlopenflags(previous_dlopenflags)
  return mod
