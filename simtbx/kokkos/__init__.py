from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
import cctbx.uctbx # possibly implicit
ext = bp.import_ext("simtbx_kokkos_ext")
from simtbx_kokkos_ext import *

