"""Automated map sharpening/blurring to optimize map"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.auto_sharpen

import sys

if (__name__ == "__main__"):

  deprecated = True
  for ok in ['--force','force']:
    if ok in sys.argv:
      sys.argv.remove(ok)
      deprecated = False
  if deprecated:
    from libtbx.utils import Sorry
    raise Sorry("phenix.auto_sharpen is replaced by phenix.map_sharpening"+
      " Add the keyword '--force' to use anyway")

  from cctbx.maptbx.auto_sharpen import run
  run(args=sys.argv[1:])

