#! /usr/bin/env python
from iotbx import reflection_file_converter
import sys
try:
  reflection_file_converter.run(sys.argv[1:])
except RuntimeError, e:
  print e
