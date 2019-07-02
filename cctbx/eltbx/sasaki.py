from __future__ import absolute_import, division, print_function
import cctbx.eltbx.fp_fdp # import dependency

import boost.python
ext = boost.python.import_ext("cctbx_eltbx_sasaki_ext")
from cctbx_eltbx_sasaki_ext import *

boost.python.inject(ext.table_iterator, boost.python.py3_make_iterator)
