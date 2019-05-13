from __future__ import division
from __future__ import print_function
import sys
import boost.python
ext = boost.python.import_ext("cctbx_sgtbx_asu_ext")
from cctbx_sgtbx_asu_ext import *

def asu_show_(asu, f=None):
  if f is None:
    f = sys.stdout
  print(asu.as_string(), file=f)


direct_space_asu.show_comprehensive_summary = asu_show_
