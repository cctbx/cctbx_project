
import sys

is_enabled = False

if (not sys.version_info[1] >= 6) :
  print "Sorry, you need at least Python 2.6 to use this module."
else :
  from libtbx.queuing_system_utils.sge.core import *
  is_enabled = True
