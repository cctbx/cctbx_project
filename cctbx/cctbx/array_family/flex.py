import scitbx.array_family.flex

import boost.python
boost.python.import_ext("cctbx_array_family_flex_ext")
from scitbx_array_family_flex_ext import *
from cctbx_array_family_flex_ext import *
import cctbx_array_family_flex_ext as ext

scitbx.array_family.flex.export_to("cctbx.array_family.flex")
