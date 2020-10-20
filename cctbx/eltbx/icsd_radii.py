from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_eltbx_icsd_radii_ext")
from cctbx_eltbx_icsd_radii_ext import *

bp.inject(ext.table_iterator, bp.py3_make_iterator)
