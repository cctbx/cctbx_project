import libtbx.forward_compatibility
import sys, os

manual_date_stamp = 20090819

def _STOP(exit_status=0):
  import sys
  f = sys._getframe(1)
  print "STOP: %s(%d)" % (f.f_code.co_filename, f.f_lineno)
  sys.exit(exit_status)
__builtins__["STOP"] = _STOP

def _numstr(values, fmt="%.6g", sep=", ", brackets=("[","]")):
  return brackets[0] + sep.join([fmt % v for v in values]) + brackets[1]
__builtins__["numstr"] = _numstr

class AutoType(object):

  def __str__(self): return "Auto"

Auto = AutoType()

class mutable(object):

  __slots__ = ["value"]

  def __init__(O, value):
    O.value = value

class slots_getstate_setstate(object):

  __slots__ = []

  def __getstate__(self):
    return dict([(name, getattr(self, name)) for name in self.__slots__])

  def __setstate__(self, state):
    for name,value in state.items(): setattr(self, name, value)

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
  if ("self" in args): del args["self"]
  else:                del args["O"]
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
    if ("self" in args): del args["self"]
    else:                del args["O"]
    del args["self"]
    for param in exclude:
      del args[param]
    self.__dict__.update(args)

class group_args(object):

  def __init__(self, **keyword_arguments):
    self.__dict__.update(keyword_arguments)

  def __call__(self):
    return self.__dict__

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


class property(object):
  """ Syntactic sugar for defining class properties for those poor souls
  who must stay compatible with older versions of Python which do not
  feature the @property decorator.

  Synopsis:

     class foo(object):

        class bar(libtbx.property):
           ''' documentation of the property
               In the following, self is the object featuring the property.
           '''
           def fget(self): # getter
           def fset(self, value): # setter
           def fdel(self): # deleter
  """

  class __metaclass__(type):

    def __new__(meta, name, bases, defs):
      if bases == (object,):
        # this is this class
        return type.__new__(meta, name, bases, defs)
      else:
        # this is some heir of this class
        return __builtins__['property'](fget=defs.get("fget"),
                                        fset=defs.get("fset"),
                                        fdel=defs.get("fdel"),
                                        doc=defs.get("__doc__"))
