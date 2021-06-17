from __future__ import absolute_import, division, print_function

import os

from libtbx.utils import escape_sh_double_quoted, gzip_open, bz2_open
from libtbx import easy_run
from libtbx.str_utils import show_string
from six.moves import cStringIO as StringIO

def for_reading(file_name, mode="r", gzip_mode="rb"):
  assert mode in ["r", "rb"]
  assert gzip_mode in ["r", "rb", "rt"]
  file_name = os.path.expanduser(file_name)
  if file_name.endswith(".gz"):
    return gzip_open(file_name=file_name, mode=gzip_mode)
  if file_name.endswith(".Z"):
    return StringIO(easy_run.fully_buffered(
      command='gunzip -c "%s"' % escape_sh_double_quoted(file_name),
      stdout_splitlines=False).raise_if_errors().stdout_buffer)
  if file_name.endswith('.bz2'):
    return bz2_open(file_name=file_name, mode=mode)
  try:
    return open(file_name, mode)
  except IOError as e:
    raise IOError(
      "Cannot open file for reading: %s\n" % show_string(file_name)
      + "  "+str(e))

def for_writing(file_name, mode="w", gzip_mode="wb"):
  assert mode in ["w", "wb", "a", "ab"]
  assert gzip_mode in ["w", "wb", "wt", "a", "ab"]
  file_name = os.path.expanduser(file_name)
  if file_name.endswith(".gz"):
    return gzip_open(file_name=file_name, mode=gzip_mode)
  elif file_name.endswith(".bz2"):
    return bz2_open(file_name=file_name, mode=mode)
  try:
    return open(file_name, mode)
  except IOError as e:
    raise IOError(
      "Cannot open file for writing: %s\n" % show_string(file_name)
      + "  "+str(e))

def file(file_name, mode):
  assert mode in ["r", "rb", "w", "wb", "a", "ab"]
  file_name = os.path.expanduser(file_name)
  if mode[0] in ("r", "rb"):
    return for_reading(file_name=file_name, mode=mode)
  return for_writing(file_name=file_name, mode=mode)

def exercise():
  import sys
  for file_name in sys.argv[1:]:
    gzip_mode = "rb"
    if file_name.endswith('.gz'):
      gzip_mode = "rt"
    assert for_reading(file_name=file_name, gzip_mode=gzip_mode).read().splitlines() \
        == ["line 1", "line 2", "the end"]
  print("line 1", file=for_writing(file_name="tmp_plain"))
  print("line 1", file=for_writing(file_name="tmp.gz", gzip_mode="wt"))
  print("OK")

if __name__ == "__main__":
  exercise()
