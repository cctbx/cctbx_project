import scitbx.stl.set # import depdendency

import boost.python
ext = boost.python.import_ext("cctbx_eltbx_chemical_elements_ext")
from cctbx_eltbx_chemical_elements_ext import *
