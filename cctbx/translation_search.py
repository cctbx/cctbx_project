import cctbx.maptbx # import dependency

import boost.python
ext = boost.python.import_ext("cctbx_translation_search_ext")
from cctbx_translation_search_ext import *
