"""Generates a multi-criterion validation kinemage file, for viewing in KiNG"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.kinemage
# LIBTBX_SET_DISPATCHER_NAME molprobity.kinemage

import sys
from mmtbx.kinemage import validation

if __name__ == "__main__":
  validation.run(sys.argv[1:])

