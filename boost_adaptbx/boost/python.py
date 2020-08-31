from __future__ import absolute_import, division, print_function

import inspect
import os
import re
import sys
import warnings

from libtbx import cpp_function_name

symbol_not_found_pat = re.compile(
  r"[Ss]ymbol[ ]not[ ]found: \s* (\w+) $", re.X | re.M | re.S)

def import_ext(name, optional=False):
  components = name.split(".")
  if (len(components) > 1):
    __import__(".".join(components[:-1]))
  previous_dlopenflags = None
  if (sys.platform.startswith("linux")) :
    previous_dlopenflags = sys.getdlopenflags()
    sys.setdlopenflags(0x100|0x2)
  try: mod = __import__(name)
  except ImportError as e:
    if (optional): return None
    error_msg = str(e)
    m = symbol_not_found_pat.search(error_msg)
    if m:
      error_msg = (  error_msg[:m.start(1)]
                   + cpp_function_name.demangle(m.group(1))
                   + error_msg[m.end(1):])
    raise ImportError(
      "\n  ".join(['__import__("%s"): %s' % (name, error_msg), "sys.path:"]
      + ["  "+p for p in sys.path]))
  for comp in components[1:]:
    mod = getattr(mod, comp)
  if (previous_dlopenflags is not None):
    sys.setdlopenflags(previous_dlopenflags)
  return mod

ext = import_ext("boost_python_meta_ext")

try: streambuf = ext.streambuf
except AttributeError: pass # XXX backward compatibility 2009-11-24
try: ostream = ext.ostream
except AttributeError: pass

if os.getenv("BOOST_ADAPTBX_ENABLE_TRACE"):
  ext.enable_signals_backtrace_if_possible()

class floating_point_exceptions_type(object):
  __shared_state = {'initialised': False}

  def __init__(self):
    self.__dict__ = self.__shared_state
    if not self.initialised:
      division_by_zero = bool(os.getenv("BOOST_ADAPTBX_TRAP_FPE", self.division_by_zero_trapped))
      invalid = bool(os.getenv("BOOST_ADAPTBX_TRAP_INVALID",self.invalid_trapped))
      overflow = bool(os.getenv("BOOST_ADAPTBX_TRAP_OVERFLOW",self.overflow_trapped))
      ext.trap_exceptions(division_by_zero, invalid, overflow)
      self.initialised = True

  def division_by_zero_trapped():
    def fget(self):
      return ext.is_division_by_zero_trapped()
    def fset(self, flag):
      if flag == self.division_by_zero_trapped: return
      ext.trap_exceptions(division_by_zero=flag,
                          invalid=self.invalid_trapped,
                          overflow=self.overflow_trapped)
    return locals()
  division_by_zero_trapped = property(**division_by_zero_trapped())

  def invalid_trapped():
    def fget(self):
      return ext.is_invalid_trapped()
    def fset(self, flag):
      if flag == self.fget(): return
      ext.trap_exceptions(division_by_zero=self.invalid_trapped,
                          invalid=flag,
                          overflow=self.overflow_trapped)
    return locals()
  invalid_trapped = property(**invalid_trapped())

  def overflow_trapped():
    def fget(self):
      return ext.is_overflow_trapped()
    def fset(self, flag):
      if flag == self.fget(): return
      ext.trap_exceptions(division_by_zero=self.overflow_trapped,
                          invalid=self.invalid_trapped,
                          overflow=flag)
    return locals()
  overflow_trapped = property(**overflow_trapped())

floating_point_exceptions = floating_point_exceptions_type()


class trapping(object):
  """ Synopsis:

      >>> import boost_adaptbx.boost.python as bp
      >>> from scitbx.array_family import flex
      >>> a = flex.double((0, 0, 0))
      >>> with bp.trapping(division_by_zero=False):
      >>>   b = 1/a
      >>> tuple(b)
      (inf, inf, inf)
      >>> 1/a
      ... CRASH ...
  """
  def __init__(self, division_by_zero=True, invalid=True, overflow=True):
    self.division_by_zero = ext.is_division_by_zero_trapped()
    self.invalid = ext.is_invalid_trapped()
    self.overflow = ext.is_overflow_trapped()
    ext.trap_exceptions(division_by_zero, invalid, overflow)


  def __enter__(self):
    return self

  def __exit__(self, exc_type, exc_val, exc_tb):
    ext.trap_exceptions(self.division_by_zero, self.invalid, self.overflow)


meta_class = ext.holder.__class__
platform_info = ext.platform_info()
assert len(platform_info) > 0 # please disable this assertion and send email to cctbx@cci.lbl.gov

def c_sizeof(typename):
  pattern = "sizeof(%s) = " % typename
  for line in platform_info.splitlines():
    if (line.startswith(pattern)):
      break
  else:
    raise RuntimeError('bp.platform_info: "%s" not found.' % pattern)
  return int(line[len(pattern):])

sizeof_void_ptr = c_sizeof("void*")


class gcc_version(object):

  def __init__(self):
    pat = r" \s* = \s* (\d+) \s+"
    m = re.search("__GNUC__ %s __GNUC_MINOR__ %s __GNUC_PATCHLEVEL__ %s"
                  % ((pat,)*3),
                  platform_info, re.X|re.M|re.S)
    if not m:
      self.major, self.minor, self.patchlevel = None, None, None
    else:
      self.major, self.minor, self.patchlevel = tuple(
        [ int(x) for x in m.groups() ])

  def __bool__(self):
    return self.major is not None

  __nonzero__ = __bool__

  def __str__(self):
    if self:
      return "%i.%i.%i" % (self.major, self.minor, self.patchlevel)
    else:
      return "GCC, it is not"


def inject(target_class, *mixin_classes):
   '''Add entries from python class dictionaries to a boost extension class.

      It is used as follows:

            class _():
              def method(...):
                ...
            bp.inject(extension_class, _)

      instead of the previous mechanism of

            class _(bp.injector, extension_class):
              def method(...):
                ...

      which does not work in python 3.
   '''
   for m in mixin_classes:
     for key, value in m.__dict__.items():
       if key not in ("__init__",
                      "__del__",
                      "__module__",
                      "__file__",
                      "__dict__") and (key != '__doc__' or value):
         setattr(target_class, key, value)

def inject_into(target_class, *mixin_classes):
  '''Add entries from python class dictionaries to a boost extension class.

     It is used as follows:

           @bp.inject_into(extension_class)
           class _():
             def method(...):
               ...
  '''
  def _inject(c):
    if inspect.isclass(c):
      inject(target_class, c, *mixin_classes)
    else:
      setattr(target_class, c.__name__, c)
      class empty_class:
        pass
      inject(target_class, empty_class, *mixin_classes)
  return _inject

def deprecate_method(boost_object, method_name):
  original_method = getattr(boost_object, method_name)

  def deprecation_helper(*args, **kwargs):
    warnings.warn(
      "The method {method_name} is deprecated and will be removed shortly".format(
        method_name=method_name
      ),
      DeprecationWarning,
      stacklevel=2,
    )
    return original_method(*args, **kwargs)

  setattr(boost_object, method_name, deprecation_helper)


def process_docstring_options(env_var="BOOST_ADAPTBX_DOCSTRING_OPTIONS"):
  from_env = os.environ.get(env_var)
  if (from_env is None): return None
  try:
    return eval(
      "docstring_options(%s)" % from_env,
      {},
      {"docstring_options": ext.docstring_options})
  except KeyboardInterrupt: raise
  except Exception as e:
    from libtbx.str_utils import show_string
    raise RuntimeError(
      "Error processing %s=%s\n" % (env_var, show_string(from_env))
      + "  Exception: %s\n" % str(e)
      + "  Valid example:\n"
      + '    %s="show_user_defined=True,show_signatures=False"' % env_var)

docstring_options = process_docstring_options()

class py3_make_iterator:
  def __next__(obj):
    return obj.next()
