# optik141 is the unmodified distribution, except for:
# mv lib optik
import sys, os
if (sys.version_info[0] < 3 or sys.version_info[1] < 3):
  sys.path.insert(0, os.path.join(os.environ["LIBTBX_DIST"], "optik141"))
  from optik import *
else:
  from optparse import *
