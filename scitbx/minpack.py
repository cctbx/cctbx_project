import scitbx.array_family.flex # import dependency

import boost.python
ext = boost.python.import_ext("scitbx_minpack_ext")
from scitbx_minpack_ext import *
