import sys
import boost.python
ext = boost.python.import_ext("cctbx_sgtbx_asu_ext")
from cctbx_sgtbx_asu_ext import *

def asu_show_(asu, f=None):
  if f is None:
    f = sys.stdout
  print >>f, asu.as_string()


direct_space_asu.show_comprehensive_summary = asu_show_
