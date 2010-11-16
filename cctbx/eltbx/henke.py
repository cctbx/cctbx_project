import cctbx.eltbx.fp_fdp # import dependency

import boost.python
ext = boost.python.import_ext("cctbx_eltbx_henke_ext")
from cctbx_eltbx_henke_ext import *
