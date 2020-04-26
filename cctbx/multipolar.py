from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex # for tuple mappings

import boost_adaptbx.python
ext = boost_adaptbx.python.import_ext("cctbx_multipolar_ext")
from cctbx_multipolar_ext import *

