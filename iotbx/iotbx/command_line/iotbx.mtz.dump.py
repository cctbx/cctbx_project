#! /bin/env python
from iotbx.mtz import dump
import sys

def run():
  for file_name in sys.argv[1:]:
    print "Reading file:", file_name
    dump.dump(file_name)

if (__name__ == "__main__"):
  run()
