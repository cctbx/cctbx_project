"""Remove outliers from data"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.remove_outliers

from mmtbx.scaling import remove_outliers
import sys

if (__name__ == "__main__"):
  remove_outliers.run(args=sys.argv[1:])

