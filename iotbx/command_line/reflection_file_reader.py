"""Read any reflection file"""
from __future__ import absolute_import, division, print_function
from iotbx import reflection_file_reader
import sys

def run():
  try:
    reflection_file_reader.run(sys.argv[1:])
  except RuntimeError as e:
    print(e)

if (__name__ == "__main__"):
  run()
