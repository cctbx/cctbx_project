from __future__ import absolute_import, division, print_function
import libtbx.forward_compatibility
import os
import sys

from libtbx.forward_compatibility import object

manual_date_stamp = 20090819

def _STOP(exit_status=0):
  f = sys._getframe(1)
  print("STOP: %s(%d)" % (f.f_code.co_filename, f.f_lineno))
  sys.exit(exit_status)
__builtins__["STOP"] = _STOP

def _numstr(
      values,
      fmt="%.6g",
      sep=", ",
      brackets=("[","]"),
      zero_threshold=None):
  flds = []
  for v in values:
    if (v is None):
      s = " "*(max(0,len(fmt % 0)-4)) + "None"
    else:
      if (zero_threshold is not None and abs(v) <= zero_threshold):
        v = 0
      s = fmt % v
      if (s.strip().replace("0", "") in ["-", "-."]):
        s = fmt % 0
    flds.append(s)
  return brackets[0] + sep.join(flds) + brackets[1]
__builtins__["numstr"] = _numstr

class AutoType(object):
  """
  Class for creating the Auto instance, which mimics the behavior of None
  with respect to the 'is' and '==' operators; this is used throughout
  CCTBX to indicate parameters that should be determined automatically.

  Examples
  --------
  >>> def f(optional=libtbx.Auto)
  ...    if optional is libtbx.Auto:
  ...        optional = 5
  ...    return optional
  ...
  >>> print(f())
  5
  >>> print(f(optional=10))
  10
  """
  singleton = None

  def __str__(self): return "Auto"
  def __eq__(self, other):
    return isinstance(other, self.__class__)
  def __ne__(self, other):
    return not self.__eq__(other)
  def __hash__(self):
    '''AutoType behaves as a singleton, so return the same hash value for all instances.'''
    return hash(AutoType)
  def __new__(cls):
    if cls.singleton is None:
      cls.singleton = super(AutoType, cls).__new__(cls)
    return cls.singleton

Auto = AutoType()

class mpi_import_guard:
  disable_mpi = False

class slots_getstate_setstate(object):
  """
  Implements getstate and setstate for classes with __slots__ defined. Allows an
  object to easily pickle only certain attributes.

  Examples
  --------
  >>> class sym_pair(libtbx.slots_getstate_setstate):
  ...     __slots__ = ["i_seq", "j_seq"]
  ...     def __init__(self, i_seq, j_seq):
  ...         self.i_seq = i_seq
  ...         self.j_seq = j_seq
  ...
  """

  __slots__ = []

  def __getstate__(self):
    """
    The name of some attributes may start with a double underscore such as
    cif_types.comp_comp_id.__rotamer_info. Python name mangling will rename such
    an attribute to _comp_comp_id_rotamer_info. Our __getstate__ function would then
    complain that the __slots__ list contains the non-existent attribute __rotamer_info.
    To fix this we manually mangle attributes with the compiler.misc.mangle function
    which does the right name mangling.
    """
    import warnings
    warning_filters = warnings.filters[:]
    show_warning = warnings.showwarning

    try:
      # avoid printing deprecation warning to stderr when loading mangle
      warnings.simplefilter("ignore")
      from libtbx.utils import mangle

    finally:
      warnings.showwarning = show_warning
      warnings.filters = warning_filters

    mnames = [ mangle(name, self.__class__.__name__) for name in self.__slots__ ]

    return dict([(name, getattr(self, name)) for name in mnames])

  def __setstate__(self, state):
    for name,value in state.items():
      if isinstance(name, bytes):
        name = name.decode('utf8')
      setattr(self, name, value)

class mutable(slots_getstate_setstate):
  __slots__ = ["value"]

  def __init__(O, value):
    O.value = value

class slots_getstate_setstate_default_initializer(slots_getstate_setstate):
  """
  Merges together functionality from slots_getstate_setstate with
  adopt_optional_init_args.

  Examples
  --------
  >>> class sym_pair(libtbx.slots_getstate_setstate_default_initializer):
  ...     __slots__ = ["i_seq", "j_seq"]
  ...
  >>> svm_pair(i_seq=1, j_seq=2)
  >>> print(svm_pair.i_seq)
  1
  """
  def __init__(self, **kwds):
    kwds = dict(kwds)
    for key in kwds :
      setattr(self, key, kwds.get(key, None))

class unpicklable(object):
  """
  An inheritable class that raises a runtime exception that an object is
  unpicklable.
  """

  def _raise_error(O):
    raise RuntimeError(
      "pickling of %s objects is disabled." % O.__class__.__name__)

  def __getinitargs__(O): O._raise_error()
  def __getstate__(O): O._raise_error()
  def __setstate__(O, state): O._raise_error()

class dict_with_default_0(dict):
  def __missing__(self, key):
    return 0

def adopt_init_args(obj, args, exclude=(), hide=False):
  """
  Adopts the initial arguments passed to an object, allowing developers to skip
  the tedious task of assigning each attribute of an instance in its __init__
  method.

  Parameters
  ----------
  obj : object
  args : list
  exclude : list of str
  hide : bool, optional

  Examples
  --------
  >>> class foo(object):
  ...     def __init__(self, x, y=1, z=None):
  ...         adopt_init_args(self, locals())
  ...
  >>> a = foo('a', z=10)
  >>> assert a.x == 'a'
  >>> assert a.y == 1
  >>> assert a.z == 10
  """
  if ("self" in args): del args["self"]
  else:                del args["O"]
  for param in exclude:
    del args[param]
  if (hide == False):
    for key in args.keys():
      if key.startswith("__") and key.endswith("__"):
        continue
      assert not hasattr(obj.__dict__, key)
    obj.__dict__.update(args)
  else:
    for key in args.keys():
      _key = "_" + key
      assert not hasattr(obj.__dict__, _key)
      obj.__dict__[_key] = args[key]

def adopt_optional_init_args(obj, kwds):
  """
  Easy management of long list of arguments with default value
  passed to __init__.

  Parameters
  ----------
  obj : object
  kwds : dict

  Examples
  --------
  >>> class foo(object):
  ...     z = 1
  ...     def __init__(self, **kwds):
  ...       libtbx.adopt_optional_init_args(self, kwds)
  ...
  >>> a = foo()
  >>> assert a.z == 1
  >>> a = foo(z=10)
  >>> assert a.z == 10
  """
  for k,v in kwds.items():
    if not hasattr(obj.__class__, k):
      raise RuntimeError("%s must be a class attribute of %s to "
                         "be adopted as optional init argument "
                         "by an instance of that class."
                         % (k, obj.__class__))
    setattr(obj, k, v)

class dda(object):
  dynamic_attributes_disabled = False
  def __setattr__(self, k, v):
    if(self.dynamic_attributes_disabled and not hasattr(self, k)):
      raise TypeError("Dynamic attributes disabled.")
    object.__setattr__(self, k, v)

  def stop_dynamic_attributes(self):
      self.dynamic_attributes_disabled = True

class group_args(dda):
  """
  Class to build an arbitrary object from a list of keyword arguments.

  Examples
  --------
  >>> from libtbx import group_args
  >>> obj = group_args(a=1, b=2, c=3)
  >>> print(obj.a, obj.b, obj.c)
  1 2 3

  Once stop_dynamic_attributes is called, adding new attributes won't be
  possible, that is this:

  obj.tmp=10

  will fail.
  """

  def __init__(self, **keyword_arguments):
    self.__dict__.update(keyword_arguments)

  def __call__(self):
    return self.__dict__

  def get(self,kw, default_value = None):
    return self.__dict__.get(kw, default_value)

  def keys(self):
    return self.__dict__.keys()

  def __repr__(self):
    outl = "group_args"
    from libtbx.utils import to_str,sys
    for attr in sorted(self.__dict__.keys()):
      tmp=getattr(self, attr)
      if (sys.version_info.major < 3) and (isinstance(tmp,unicode)):
        tmp = to_str(tmp)
      if str(tmp).find("ext.atom ")>-1:
        outl += "\n  %-30s : %s" % (attr, tmp.quote())
      else:
        outl += "\n  %-30s : %s" % (attr, tmp)
    return outl

  def merge(self, other):
    """ To merge other group_args into self.
    Overwrites matching fields!!!"""
    self.__dict__.update(other.__dict__)

  def add_if_missing(self, other, add_if_self_is_none = False):
    """ takes values from other only if not present at all in self
      Optionally add if value in self is None"""
    self_keys = list(self.keys())
    for key in other.keys():
      if key.startswith("__"): continue
      if (not (key in self_keys)) or (add_if_self_is_none and
          (self.get(key) is None)):
        self.add(key,other.get(key))

  def add(self,key=None,value=None):
    self.__dict__[key]=value

  def delete(self, key = None):
    if key in self.keys():
      self.__dict__[key] = None
  def copy(self):
    """ produce shallow copy of self by converting to dict and back"""
    return group_args(**self().copy())

if os.environ.get("LIBTBX_PRINT_TRACE"):
  import libtbx.start_print_trace

if sys.platform == "cygwin":
  # work around cygwin problem: open() doesn't work on symbolic links
  builtin_open = __builtins__["open"]
  def open_realpath(name, mode="r", buffering=-1):
    try: return builtin_open(name, mode, buffering)
    except KeyboardInterrupt: raise
    except Exception: pass
    name = os.path.realpath(name)
    return builtin_open(name, mode, buffering)
  __builtins__["open"] = open_realpath
  __builtins__["file"] = open_realpath
