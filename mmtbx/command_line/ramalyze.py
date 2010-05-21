# LIBTBX_SET_DISPATCHER_NAME phenix.ramalyze

import sys
from mmtbx.validation.ramalyze import ramalyze

if __name__ == "__main__":
  r = ramalyze()
  output_list, coot_todo_list = r.run(sys.argv[1:])
