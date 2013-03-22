from __future__ import division
from scitbx.array_family import flex # import dependency
try:
  import boost.python
except Exception:
  ext = None
else:
  ext = boost.python.import_ext("dxtbx_ext", optional = True)

if not ext is None:
  from dxtbx_ext import *

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

