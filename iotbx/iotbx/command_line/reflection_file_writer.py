#! /usr/bin/env python
from iotbx import reflection_file_writer
import sys
try:
  reflection_file_writer.run(sys.argv[1:])
except RuntimeError, e:
  print e
