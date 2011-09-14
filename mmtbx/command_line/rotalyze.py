# LIBTBX_SET_DISPATCHER_NAME phenix.rotalyze

import sys
from mmtbx.validation.rotalyze import rotalyze

if __name__ == "__main__":
  r = rotalyze()
  output_list, coot_todo_list = r.run(sys.argv[1:])
