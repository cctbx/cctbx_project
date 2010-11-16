import cctbx.eltbx.fp_fdp # import dependency

import boost.python
ext = boost.python.import_ext("cctbx_eltbx_sasaki_ext")
from cctbx_eltbx_sasaki_ext import *
