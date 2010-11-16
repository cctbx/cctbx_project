from scitbx.array_family import flex # import dependency

import boost.python
ext = boost.python.import_ext("scitbx_r3_utils_ext")
from scitbx_r3_utils_ext import *
