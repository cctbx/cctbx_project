#! /usr/bin/env python
from iotbx import reflection_file_reader
import sys
if (len(sys.argv) == 1):
  print reflection_file_reader.usage()
else:
  try:
    reflection_file_reader.run(sys.argv[1:])
  except RuntimeError, e:
    print e
