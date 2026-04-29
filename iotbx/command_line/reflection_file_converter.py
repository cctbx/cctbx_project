"""Convert or modify reflection file"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.reflection_file_converter
# LIBTBX_SET_DISPATCHER_NAME iotbx.reflection_file_converter

from iotbx import reflection_file_converter
import sys

def run():
  try:
    reflection_file_converter.run(args=sys.argv[1:])
  except RuntimeError as e:
    print(e)

if (__name__ == "__main__"):
  run()
