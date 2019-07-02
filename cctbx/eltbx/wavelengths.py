from __future__ import absolute_import, division, print_function
import boost.python
ext = boost.python.import_ext("cctbx_eltbx_wavelengths_ext")
from cctbx_eltbx_wavelengths_ext import *

boost.python.inject(ext.characteristic_iterator, boost.python.py3_make_iterator)
