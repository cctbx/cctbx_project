import cctbx.array_family.flex

import boost.python
ext = boost.python.import_ext("iotbx_mtz_wrapper_ext")
from iotbx_mtz_wrapper_ext import *
import iotbx_mtz_wrapper_ext as ext
