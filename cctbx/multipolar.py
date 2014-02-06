from __future__ import division
from cctbx.array_family import flex # for tuple mappings

import boost.python
ext = boost.python.import_ext("cctbx_multipolar_ext")
from cctbx_multipolar_ext import *

