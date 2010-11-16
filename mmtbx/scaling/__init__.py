import cctbx.array_family.flex # import dependency

import boost.python
ext = boost.python.import_ext("mmtbx_scaling_ext")
from mmtbx_scaling_ext import *
