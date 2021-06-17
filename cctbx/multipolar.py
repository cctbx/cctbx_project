from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex # for tuple mappings

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_multipolar_ext")
from cctbx_multipolar_ext import *
