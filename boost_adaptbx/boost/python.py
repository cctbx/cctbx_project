import sys, os
import re
from libtbx import cpp_function_name
symbol_not_found_pat = re.compile(
  r"[Ss]ymbol[ ]not[ ]found: \s* (\w+) $", re.X | re.M | re.S)

python_libstdcxx_so = None
if (sys.platform.startswith("linux")):
  from libtbx import easy_run
  for line in easy_run.fully_buffered(
                command='/usr/bin/ldd "%s"' % sys.executable).stdout_lines:
    if (line.strip().startswith("libstdc++.so")):
      python_libstdcxx_so = line.split()[0]
      break

def import_ext(name, optional=False):
  components = name.split(".")
  if (len(components) > 1):
    __import__(".".join(components[:-1]))
  previous_dlopenflags = None
  if (sys.platform.startswith("linux")) :
    previous_dlopenflags = sys.getdlopenflags()
    sys.setdlopenflags(0x100|0x2)
  try: mod = __import__(name)
  except ImportError, e:
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
  if (python_libstdcxx_so is not None):
    mod_file = getattr(mod, "__file__", None)
    if (mod_file is not None):
      for line in easy_run.fully_buffered(
                    command='/usr/bin/ldd "%s"' % mod_file).stdout_lines:
        if (line.strip().startswith("libstdc++.so")):
          mod_libstdcxx_so = line.split()[0]
          if (mod_libstdcxx_so != python_libstdcxx_so):
            raise SystemError("""\
FATAL: libstdc++.so mismatch:
  %s: %s
  %s: %s""" % (sys.executable, python_libstdcxx_so,
               mod_file, mod_libstdcxx_so))
          break
  return mod

ext = import_ext("boost_python_meta_ext")

try: streambuf = ext.streambuf
except AttributeError: pass # XXX backward compatibility 2009-11-24
try: ostream = ext.ostream
except AttributeError: pass

if ("BOOST_ADAPTBX_SIGNALS_DEFAULT" not in os.environ):
  ext.enable_signals_backtrace_if_possible()


class floating_point_exceptions_type(object):

  __shared_state = {'initialised': False}

  def __init__(self, division_by_zero, invalid, overflow):
    self.__dict__ = self.__shared_state
    if not self.initialised:
      if "BOOST_ADAPTBX_FPE_DEFAULT" in os.environ:
        division_by_zero = self.division_by_zero_trapped
        invalid = self.invalid_trapped
        overflow = self.overflow_trapped
      elif "BOOST_ADAPTBX_FE_DIVBYZERO_DEFAULT" in os.environ:
        division_by_zero = self.division_by_zero_trapped
      elif "BOOST_ADAPTBX_FE_INVALID_DEFAULT" in os.environ:
        invalid = self.invalid_trapped
      elif "BOOST_ADAPTBX_FE_OVERFLOW_DEFAULT" in os.environ:
        overflow = self.overflow_trapped
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

def floating_point_exceptions():
  import libtbx.load_env
  if (libtbx.env.is_development_environment()):
    flag = True
  else:
    flag = False
  return floating_point_exceptions_type(
    division_by_zero=flag, invalid=flag, overflow=flag)
floating_point_exceptions = floating_point_exceptions()

meta_class = ext.holder.__class__
platform_info = ext.platform_info()
assert len(platform_info) > 0 # please disable this assertion and send email to cctbx@cci.lbl.gov

def c_sizeof(typename):
  pattern = "sizeof(%s) = " % typename
  for line in platform_info.splitlines():
    if (line.startswith(pattern)):
      break
  else:
    raise RuntimeError('boost.python.platform_info: "%s" not found.' % pattern)
  return int(line[len(pattern):])

sizeof_void_ptr = c_sizeof("void*")


class gcc_version(object):

  def __init__(self):
    import re
    pat = r" \s* = \s* (\d+) \s+"
    m = re.search("__GNUC__ %s __GNUC_MINOR__ %s __GNUC_PATCHLEVEL__ %s"
                  % ((pat,)*3),
                  platform_info, re.X|re.M|re.S)
    if not m:
      self.major, self.minor, self.patchlevel = None, None, None
    else:
      self.major, self.minor, self.patchlevel = tuple(
        [ int(x) for x in m.groups() ])

  def __nonzero__(self):
    return self.major is not None

  def __str__(self):
    if self:
      return "%i.%i.%i" % (self.major, self.minor, self.patchlevel)
    else:
      return "GCC, it is not"


class injector(object):

  class __metaclass__(meta_class):

    def __init__(self, name, bases, dict):
      if (len(bases) > 1):
        # bases[0] is this injector
        target = bases[1] # usually a Boost.Python class
        def setattr_from_dict(d):
          for k,v in d.items():
            if (k not in (
                  "__init__",
                  "__del__",
                  "__module__",
                  "__file__",
                  "__dict__")):
              if (k != "__doc__" or v is not None):
                setattr(target,k,v)
        setattr_from_dict(dict)
        for b in bases[2:]: # usually mix-in classes, if any
          setattr_from_dict(b.__dict__)
      return type.__init__(self, name, (), {})

def process_docstring_options(env_var="BOOST_ADAPTBX_DOCSTRING_OPTIONS"):
  from_env = os.environ.get(env_var)
  if (from_env is None): return None
  try:
    return eval(
      "docstring_options(%s)" % from_env,
      {},
      {"docstring_options": ext.docstring_options})
  except KeyboardInterrupt: raise
  except Exception, e:
    from libtbx.str_utils import show_string
    raise RuntimeError(
      "Error processing %s=%s\n" % (env_var, show_string(from_env))
      + "  Exception: %s\n" % str(e)
      + "  Valid example:\n"
      + '    %s="show_user_defined=True,show_signatures=False"' % env_var)

docstring_options = process_docstring_options()
