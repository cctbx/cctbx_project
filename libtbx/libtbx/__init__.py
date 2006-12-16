import libtbx.forward_compatibility
import sys, os

if (getattr(sys, "api_version", 0) >= 1013):

  class dict_with_default_0(dict):

    def __missing__(self, key):
      return 0

else:

  class dict_with_default_0(dict):

    def __getitem__(self, key):
      try: return dict.__getitem__(self, key)
      except KeyError: pass
      val = 0
      dict.__setitem__(self, key, val)
      return val

def adopt_init_args(obj, args, exclude=(), hide=False):
  del args["self"]
  for param in exclude:
    del args[param]
  if (hide == False):
    for key in args.keys():
      assert not hasattr(obj.__dict__, key)
    obj.__dict__.update(args)
  else:
    for key in args.keys():
      _key = "_" + key
      assert not hasattr(obj.__dict__, _key)
      obj.__dict__[_key] = args[key]

if (os.environ.has_key("LIBTBX_PRINT_TRACE")):
  import libtbx.start_print_trace
