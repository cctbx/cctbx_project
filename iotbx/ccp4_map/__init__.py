import cctbx.array_family.flex # import dependency

import boost.python
ext = boost.python.import_ext("iotbx_ccp4_map_ext")
from iotbx_ccp4_map_ext import *
import iotbx_ccp4_map_ext as ext
