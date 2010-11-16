from cctbx.array_family import flex # import dependency

import boost.python
ext = boost.python.import_ext("rstbx_ext")
from rstbx_ext import *
