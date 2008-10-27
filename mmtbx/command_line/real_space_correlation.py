from mmtbx import real_space_correlation
import sys

if(__name__ == "__main__"):
  real_space_correlation.cmd_run(args         = sys.argv[1:],
                                 command_name = "mmtbx.real_space_correlation")
