from __future__ import division
# Not sure if the next line is needed (copied this script from mmtbx/tls/__init__.py)
#import cctbx.array_family.flex # import dependency

import boost.python
ext = boost.python.import_ext("mmtbx_tls_decompose_ext")
from mmtbx_tls_decompose_ext import *
