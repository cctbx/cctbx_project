import scitbx.array_family.flex

import boost.python
ext = boost.python.import_ext("cctbx_array_family_flex_ext")
from scitbx_array_family_flex_ext import *
from cctbx_array_family_flex_ext import *

scitbx.array_family.flex.export_to("cctbx.array_family.flex")
