import cctbx.array_family.flex # import dependency

import boost.python
ext = boost.python.import_ext("mmtbx_tls_ext")
from mmtbx_tls_ext import *
