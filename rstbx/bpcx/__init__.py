from __future__ import absolute_import, division, print_function
# usage:
#
# from rstbx.bpcx import sensor

from cctbx.array_family import flex # import dependency
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("rstbx_bpcx_detector_model_ext")
from rstbx_bpcx_detector_model_ext import *
