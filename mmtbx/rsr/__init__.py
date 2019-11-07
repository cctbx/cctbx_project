from __future__ import absolute_import, division, print_function
import cctbx.array_family.flex # import dependency

import boost.python
ext = boost.python.import_ext("mmtbx_rsr_ext")
from mmtbx_rsr_ext import *
