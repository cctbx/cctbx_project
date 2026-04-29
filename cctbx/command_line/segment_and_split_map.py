"""Segment a map and write out one map per segment"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.segment_and_split_map

import sys

if (__name__ == "__main__"):
  from cctbx.maptbx.segment_and_split_map import run
  run(args=sys.argv[1:])

