from __future__ import absolute_import, division, print_function

import logging

try:
  import boost.python
except Exception:
  ext = None
else:
  ext = boost.python.import_ext("dxtbx_ext", optional = True)

if ext is not None:
  from dxtbx_ext import *

logging.getLogger('dxtbx').addHandler(logging.NullHandler())

class IncorrectFormatError(RuntimeError):
  '''
  An exception class for an incorrect format error
  '''
  def __init__(self, format_instance, filename):
    super(IncorrectFormatError, self).__init__(
      "Could not open %s as %s" % (filename, str(format_instance)))


def load(filename):
  """Use DXTBX to load the input filename.

  :param filename:  The input filename
  :type  filename:  os.PathLike or str or bytes
  :returns:         A dxtbx Format-subclass instance for the file type
  :raises IOError:  if the file format could not be determined
  """
  from dxtbx.format.Registry import Registry
  # Unwrap PEP-519-style objects. This covers py.path, pathlib, ...
  if hasattr(filename, "__fspath__"):
    filename = filename.__fspath__()
  format_instance = Registry.find(filename)
  return format_instance(filename)

def get_print():
  '''Get the builtin print function.'''
  try:
    import __builtin__
  except ImportError:
    import builtins as __builtin__
  return __builtin__.print
