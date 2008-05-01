import libtbx.forward_compatibility
import sys, os

class AutoType(object):

  def __str__(self): return "Auto"

Auto = AutoType()

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

def adopt_optional_init_args(obj, kwds):
  """\
  Description:
    Easy management of long list of arguments with default value
    passed to __init__
  Synopsis:
    class foo(object):
      z = 1
      def __init__(self, x, y, **kwds):
        self.x = x
        self.y = y
        adopt_optional_init_args(self, kwds)

    a = foo(x,y)
    assert a.z == 1
    a = foo(x,y, z=10)
    assert a.z == 10
  """
  for k,v in kwds.iteritems():
    if not hasattr(obj.__class__, k):
      raise RuntimeError("%s must be a class attribute of %s to "
                         "be adopted as optional init argument "
                         "by an instance of that class."
                         % (k, obj.__class__))
    setattr(obj, k, v)

class copy_init_args(object):

  def __init__(self, args, exclude=()):
    del args["self"]
    for param in exclude:
      del args[param]
    self.__dict__.update(args)

class group_args(object):

  def __init__(self, **keyword_arguments):
    self.__dict__.update(keyword_arguments)

if (os.environ.has_key("LIBTBX_PRINT_TRACE")):
  import libtbx.start_print_trace

if (sys.platform == "cygwin"):
  # work around cygwin problem: open() doesn't work on symbolic links
  builtin_open = __builtins__["open"]
  def open_realpath(name, mode="r", buffering=-1):
    try: return builtin_open(name, mode, buffering)
    except KeyboardInterrupt: raise
    except: pass
    name = os.path.realpath(name)
    return builtin_open(name, mode, buffering)
  __builtins__["open"] = open_realpath
  __builtins__["file"] = open_realpath
