from mmtbx.scaling import remove_outliers 
import sys

if (__name__ == "__main__"):
  remove_outliers.run(command_name="mmtbx.remove_outliers", args=sys.argv[1:])
