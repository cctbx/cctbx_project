from __future__ import absolute_import, division, print_function
import scitbx.stl.vector # import dependency

import boost_adaptbx.python
ext = boost_adaptbx.python.import_ext("scitbx_stl_map_ext")
from scitbx_stl_map_ext import *
