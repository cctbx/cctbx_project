#! /usr/bin/env python
from iotbx import reflection_file_reader
import sys
try:
  reflection_file_reader.run(sys.argv[1:])
except RuntimeError, e:
  print e
