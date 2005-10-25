from mmtbx.scaling import xtriage
import sys

if (__name__ == "__main__"):
  xtriage.run(command_name="mmtbx.xtriage", args=sys.argv[1:])
