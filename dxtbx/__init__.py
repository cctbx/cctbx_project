from __future__ import division
from scitbx.array_family import flex # import dependency
try:
  import boost.python
except Exception:
  ext = None
else:
  ext = boost.python.import_ext("dxtbx_ext", optional = False)

if not ext is None:
  from dxtbx_ext import *

