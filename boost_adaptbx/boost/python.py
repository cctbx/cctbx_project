import sys

python_libstdcxx_so = None
if (sys.platform.startswith("linux")):
  from libtbx import easy_run
  for line in easy_run.fully_buffered(
                command='/usr/bin/ldd "%s"' % sys.executable).stdout_lines:
    if (line.strip().startswith("libstdc++.so")):
      python_libstdcxx_so = line.split()[0]
      break

def import_ext(name):
  components = name.split(".")
  if (len(components) > 1):
    __import__(".".join(components[:-1]))
  previous_dlopenflags = None
  if (sys.platform == "linux2"):
    previous_dlopenflags = sys.getdlopenflags()
    sys.setdlopenflags(0x100|0x2)
  try: mod = __import__(name)
  except ImportError, e:
    raise ImportError(
      "\n  ".join(['__import__("%s"): %s' % (name, str(e)), "sys.path:"]
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
meta_class = ext.holder.__class__
platform_info = ext.platform_info()
assert len(platform_info) > 0 # please disable this assertion and send email to cctbx@cci.lbl.gov
sizeof_void_ptr = ext.sizeof_void_ptr()

class injector(object):
  "see boost/libs/python/doc/tutorial/doc/quickstart.txt"

  class __metaclass__(meta_class):

    def __init__(self, name, bases, dict):
      for b in bases:
        if type(b) not in (self, type):
          for k,v in dict.items():
            if (k in ("__init__", "__del__", "__module__", "__file__")):
              continue
            setattr(b,k,v)
      return type.__init__(self, name, bases, dict)
