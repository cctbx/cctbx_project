from __future__ import absolute_import, division, print_function
import sys
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_sgtbx_asu_ext")
from cctbx_sgtbx_asu_ext import *

def asu_show_(asu, f=None):
  if f is None:
    f = sys.stdout
  print(asu.as_string(), file=f)


direct_space_asu.show_comprehensive_summary = asu_show_
