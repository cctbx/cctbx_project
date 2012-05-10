# usage:
#
# from rstbx.bpcx import sensor

import boost.python
ext = boost.python.import_ext("rstbx_bpcx_detector_model_ext")
from rstbx_bpcx_detector_model_ext import *
