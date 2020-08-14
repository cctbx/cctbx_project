from __future__ import absolute_import, division, print_function
import scitbx.stl.vector # import dependency

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("scitbx_stl_map_ext")
from scitbx_stl_map_ext import *
