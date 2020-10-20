from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_eltbx_wavelengths_ext")
from cctbx_eltbx_wavelengths_ext import *

bp.inject(ext.characteristic_iterator, bp.py3_make_iterator)
