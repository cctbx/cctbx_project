from iotbx.scalepack import reader
import sys

def run():
  reader.run(sys.argv[1:])

if (__name__ == "__main__"):
  run()
