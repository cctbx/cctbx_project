import scitbx.array_family.flex
import cctbx.array_family.flex

import boost.python
boost.python.import_ext("cctbx_array_family_flex_ext")
from scitbx_array_family_flex_ext import *
from cctbx_array_family_flex_ext import *
from rstbx_array_family_flex_ext import *
import rstbx_array_family_flex_ext as ext # implicit import

scitbx.array_family.flex.export_to("rstbx.array_family.flex")
