from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.cablam_validate

import sys
try:
  from mmtbx.cablam import cablam_validate
except Exception as e:
  print(e)

if __name__ == "__main__":
  cablam_validate.run(sys.argv[1:])
