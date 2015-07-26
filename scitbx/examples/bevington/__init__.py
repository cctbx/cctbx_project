from __future__ import division
from scitbx.lstbx import normal_eqns # import dependency
import boost.python
ext = boost.python.import_ext("scitbx_examples_bevington_ext")
from scitbx_examples_bevington_ext import *
