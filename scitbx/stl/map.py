import scitbx.stl.vector # import dependency

import boost.python
ext = boost.python.import_ext("scitbx_stl_map_ext")
from scitbx_stl_map_ext import *
