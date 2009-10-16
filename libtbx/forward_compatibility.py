from __future__ import generators
import sys, os
if not hasattr(os.path, "devnull"):
  # Python 2.3 compatibility
  if os.name == "nt":
    os.path.devnull = "nul"
  else:
    os.path.devnull = "/dev/null"

if ("sorted" not in __builtins__):
  # Python 2.3 compatibility
  def sorted(iterable, cmp=None, key=None, reverse=False):
    """\
sorted(iterable, cmp=None, key=None, reverse=False) --> new sorted list"""
    assert key is None, "Not implemented."
    result = list(iterable)
    if (cmp is None): result.sort()
    else:             result.sort(cmp)
    if (reverse): result.reverse()
    return result
  __builtins__["sorted"] = sorted

if ("reversed" not in __builtins__):
  # Python 2.3 compatibility
  def reversed(seq):
    """ Return a reverse iterator. seq must be an object which supports
    the sequence protocol (the __len__() method and the __getitem__() method
    with integer arguments starting at 0). New in version 2.4. """
    i = len(seq)
    while i > 0:
      i -= 1
      yield seq[i]
  __builtins__["reversed"] = reversed

vers_info = sys.version_info[:2]
if (vers_info == (2,3) and "set" not in __builtins__):
  # Python 2.3 compatibility
  import sets
  __builtins__["set"] = sets.Set
  __builtins__["frozenset"] = sets.ImmutableSet
if (vers_info in [(2,3),(2,4),(2,5)] and not hasattr(frozenset, "isdisjoint")):
  # Python 2.3, 2.4, 2.5 compatibility
  class forward_compatibility_set_mixin(object):
    def isdisjoint(self, other):
      for value in other:
        if value in self:
          return False
      return True
  class forward_compatibility_frozenset(
        forward_compatibility_set_mixin, frozenset): pass
  class forward_compatibility_set(
        forward_compatibility_set_mixin, set): pass
  __builtins__["frozenset"] = forward_compatibility_frozenset
  __builtins__["set"] = forward_compatibility_set

class _advertise_subprocess(object):

  def __init__(self, function_id, target):
    self.function_id = function_id
    self.target = target

  def __call__(self, *args, **kwargs):
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

_install_advertise_subprocess()
