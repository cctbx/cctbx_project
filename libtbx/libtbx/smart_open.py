from libtbx.utils import escape_sh_double_quoted
from libtbx import easy_run
from libtbx.str_utils import show_string
try: import gzip
except ImportError: gzip = None
from cStringIO import StringIO

def for_reading(file_name, mode="r", gzip_mode="rb"):
  assert mode in ["r", "rb"]
  assert gzip_mode in ["r", "rb"]
  if (file_name.endswith(".gz")):
    if (gzip is None):
      raise RuntimeError(
        "gzip module not available: cannot uncompress file %s"
          % show_string(file_name))
    return gzip.open(file_name, gzip_mode)
  if (file_name.endswith(".Z")):
    return StringIO(easy_run.fully_buffered(
      command='gunzip -c "%s"' % escape_sh_double_quoted(file_name),
      stdout_splitlines=False).raise_if_errors().stdout_buffer)
  try:
    return open(file_name, mode)
  except IOError, e:
    raise IOError(
      "Cannot open file for reading: %s\n" % show_string(file_name)
      + "  "+str(e))

def exercise():
  import sys
  for file_name in sys.argv[1:]:
    assert for_reading(file_name=file_name).read().splitlines() \
        == ["line 1", "line 2", "the end"]
  print "OK"

if (__name__ == "__main__"):
  exercise()
