import scitbx.array_family.flex # import dependency

import boost.python
ext = boost.python.import_ext("scitbx_fftpack_ext")
from scitbx_fftpack_ext import *
