from scitbx.array_family import flex

import libtbx.boost_python
ext = libtbx.boost_python.import_ext("iotbx_detectors_ext")
from iotbx_detectors_ext import *
