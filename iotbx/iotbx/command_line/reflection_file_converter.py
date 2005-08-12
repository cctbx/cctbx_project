from iotbx import reflection_file_converter
import sys

def run():
  try:
    reflection_file_converter.run(sys.argv[1:])
  except RuntimeError, e:
    print e

if (__name__ == "__main__"):
  run()
