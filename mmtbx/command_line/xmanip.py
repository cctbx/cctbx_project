"""Experimental tool for manipulation of reflection data"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.xmanip

from mmtbx import xmanip
import sys

if (__name__ == "__main__"):
  xmanip.run(args=sys.argv[1:])

