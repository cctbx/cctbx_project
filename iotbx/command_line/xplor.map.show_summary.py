import iotbx.xplor.map
import sys

def run(args):
  for file_name in args:
    print "file name:", file_name
    iotbx.xplor.map.reader(file_name=file_name).show_summary(prefix="  ")
    print

if (__name__ == "__main__"):
  run(sys.argv[1:])
