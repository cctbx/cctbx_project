from __future__ import absolute_import, division
from scitbx.array_family import flex # import dependency
try:
  import boost.python
except Exception:
  ext = None
else:
  ext = boost.python.import_ext("dxtbx_ext", optional = True)

if ext is not None:
  from dxtbx_ext import *

import dxtbx.imageset # implicit import

class IncorrectFormatError(RuntimeError):
  '''
  An exception class for an incorrect format error

  '''

  def __init__(self, format_instance, filename):
    super(IncorrectFormatError, self).__init__(
      "Could not open %s as %s" % (filename, str(format_instance)))


def load(filename):
  """Use DXTBX to get the files from the input filename.

  Params:
      filename The input filename

  Returns:
      The dxtbx format instance

  """
  from dxtbx.format.Registry import Registry
  format_instance = Registry.find(filename)
  return format_instance(filename)
