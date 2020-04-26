from __future__ import absolute_import, division, print_function
from scitbx.examples import bevington # import dependency
import boost_adaptbx.python
ext = boost_adaptbx.python.import_ext("cctbx_large_scale_merging_ext")
from cctbx_large_scale_merging_ext import *
