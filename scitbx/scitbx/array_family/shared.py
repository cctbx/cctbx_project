import scitbx.array_family.flex # for slice converters

import boost.python
boost.python.import_ext("scitbx_array_family_shared_ext")
from scitbx_array_family_shared_ext import *
import scitbx_array_family_shared_ext as ext
