from __future__ import absolute_import, division, print_function
import sys
import os

_vers_info = sys.version_info[:2]

class newobject(object):
  """
  A magical object class that provides Python 2 compatibility methods::
      next
      __unicode__
      __nonzero__

  Subclasses of this class can merely define the Python 3 methods (__next__,
  __str__, and __bool__).
  """
  # verbatim copy of builtin.object from future package + __slots__ fix
  # cf https://github.com/PythonCharmers/python-future
  #    https://github.com/PythonCharmers/python-future/pull/337

  def next(self):
      if hasattr(self, '__next__'):
          return type(self).__next__(self)
      raise TypeError('newobject is not an iterator')

  def __unicode__(self):
      # All subclasses of the builtin object should have __str__ defined.
      # Note that old-style classes do not have __str__ defined.
      if hasattr(self, '__str__'):
          s = type(self).__str__(self)
      else:
          s = str(self)
      if isinstance(s, unicode):
          return s
      else:
          return s.decode('utf-8')

  def __nonzero__(self):
      if hasattr(self, '__bool__'):
          return type(self).__bool__(self)
      if hasattr(self, '__len__'):
          return type(self).__len__(self)
      # object has no __nonzero__ method
      return True

  # Are these ever needed?
  # def __div__(self):
  #     return self.__truediv__()

  # def __idiv__(self, other):
  #     return self.__itruediv__(other)

  def __long__(self):
      if not hasattr(self, '__int__'):
          return NotImplemented
      return self.__int__()  # not type(self).__int__(self)

  # def __new__(cls, *args, **kwargs):
  #     """
  #     dict() -> new empty dictionary
  #     dict(mapping) -> new dictionary initialized from a mapping object's
  #         (key, value) pairs
  #     dict(iterable) -> new dictionary initialized as if via:
  #         d = {}
  #         for k, v in iterable:
  #             d[k] = v
  #     dict(**kwargs) -> new dictionary initialized with the name=value pairs
  #         in the keyword argument list.  For example:  dict(one=1, two=2)
  #     """

  #     if len(args) == 0:
  #         return super(newdict, cls).__new__(cls)
  #     elif type(args[0]) == newdict:
  #         return args[0]
  #     else:
  #         value = args[0]
  #     return super(newdict, cls).__new__(cls, value)

  def __native__(self):
      """
      Hook for the future.utils.native() function
      """
      return object(self)

  __slots__ = []

if _vers_info < (3,0):
  # export newobject as object
  object = newobject
  del newobject
  try:
    # If builtins exists, monkeypatch it.
    import builtins
    builtins.object = object
  except ImportError:
    # builtins may be missing, eg. in case of running an installer script
    pass
else:
  object = object # Statement required to export the name

class _advertise_subprocess(object):

  def __init__(self, function_id, target):
    self.function_id = function_id
    self.target = target

  def __call__(self, *args, **kwargs):
    if ("LIBTBX_BUILD" in os.environ):
      import libtbx.load_env
      if (   libtbx.env.full_testing
          or libtbx.env.is_development_environment()):
        def is_exception():
          from sys import _getframe
          frames_back = 1
          while True:
            try: f = _getframe(frames_back)
            except ValueError: break
            else:
              p = f.f_code.co_filename
              while True:
                d,b = os.path.split(p)
                if (len(d) == 0 or len(b) == 0): break
                b = b.lower()
                if (   b.startswith("python")
                    or b == "lib"):
                  if (not os.path.isdir(p)): # non-existing path names in .pyc
                    return True              # simply give up to avoid noise
                  for n in ["os","StringIO","UserDict","UserList","UserString"]:
                    if (not os.path.isfile(os.path.join(p, n+".py"))):
                      break
                  else:
                    return True
                elif (b == "scons"):
                  if (os.path.isfile(os.path.join(p, "SConsign.py"))):
                    return True
                p = d
              frames_back += 1
          return False
        if (not is_exception()):
          from warnings import warn
          warn(
            message="%s is not safe: please use the subprocess module"
                    " or libtbx.easy_run instead." % self.function_id,
            stacklevel=2)
    return self.target(*args, **kwargs)

def _install_advertise_subprocess():
  for fn in ["system",
             "popen", "popen2", "popen3", "popen4",
             "spawnl", "spawnle", "spawnlp", "spawnlpe",
             "spawnv", "spawnve", "spawnvp", "spawnvpe"]:
    f = getattr(os, fn, None)
    if (f is not None):
      w = _advertise_subprocess(function_id="os."+fn+"()", target=f)
      setattr(os, fn, w)

if _vers_info < (3,0):
  _install_advertise_subprocess()
