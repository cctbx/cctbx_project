from __future__ import absolute_import, division, print_function
import boost_adaptbx.python
ext = boost_adaptbx.python.import_ext("cctbx_eltbx_wavelengths_ext")
from cctbx_eltbx_wavelengths_ext import *

boost_adaptbx.python.inject(ext.characteristic_iterator, boost_adaptbx.python.py3_make_iterator)
