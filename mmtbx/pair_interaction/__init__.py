from __future__ import division
import cctbx.array_family.flex # import dependency

import boost.python
ext = boost.python.import_ext("mmtbx_interaction_ext")
from mmtbx_interaction_ext import *
