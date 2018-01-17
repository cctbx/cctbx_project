from __future__ import absolute_import, division
import sys, os

vers_info = sys.version_info[:2]

if (vers_info < (2,6)):
  import cmath
  import math
  cmath.phase = lambda z: math.atan2(z.imag, z.real)
  cmath.polar = lambda z: (abs(z), cmath.phase(z))
  cmath.rect = lambda r, phi: (r*math.cos(phi), r*math.sin(phi))

  import warnings
  # This is copied verbatim from Python 2.7 library
  class catch_warnings(object):

      """A context manager that copies and restores the warnings filter upon
      exiting the context.

      The 'record' argument specifies whether warnings should be captured by a
      custom implementation of warnings.showwarning() and be appended to a list
      returned by the context manager. Otherwise None is returned by the context
      manager. The objects appended to the list are arguments whose attributes
      mirror the arguments to showwarning().

      The 'module' argument is to specify an alternative module to the module
      named 'warnings' and imported under that name. This argument is only useful
      when testing the warnings module itself.

      """

      def __init__(self, record=False, module=None):
          """Specify whether to record warnings and if an alternative module
          should be used other than sys.modules['warnings'].

          For compatibility with Python 3.0, please consider all arguments to be
          keyword-only.

          """
          self._record = record
          if (module is None) :
            self._module = sys.modules['warnings']
          else :
            self._module = module
          self._entered = False

      def __repr__(self):
          args = []
          if self._record:
              args.append("record=True")
          if self._module is not sys.modules['warnings']:
              args.append("module=%r" % self._module)
          name = type(self).__name__
          return "%s(%s)" % (name, ", ".join(args))

      def __enter__(self):
          if self._entered:
              raise RuntimeError("Cannot enter %r twice" % self)
          self._entered = True
          self._filters = self._module.filters
          self._module.filters = self._filters[:]
          self._showwarning = self._module.showwarning
          if self._record:
              log = []
              def showwarning(*args, **kwargs):
                  log.append(WarningMessage(*args, **kwargs))
              self._module.showwarning = showwarning
              return log
          else:
              return None

      def __exit__(self, *exc_info):
          if not self._entered:
              raise RuntimeError("Cannot exit %r without entering first" % self)
          self._module.filters = self._filters
          self._module.showwarning = self._showwarning
  # end of the copy
  warnings.catch_warnings = catch_warnings

  if (vers_info in [(2,3),(2,4),(2,5)]
        and not hasattr(frozenset, "isdisjoint")):
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

del vers_info

class _advertise_subprocess(object):

  def __init__(self, function_id, target):
    self.function_id = function_id
    self.target = target

  def __call__(self, *args, **kwargs):
    if (os.environ.has_key("LIBTBX_BUILD")):
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
