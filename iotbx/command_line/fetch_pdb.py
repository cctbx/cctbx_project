# LIBTBX_SET_DISPATCHER_NAME phenix.fetch_pdb

import sys

if __name__ == "__main__" :
  import iotbx.pdb.fetch
  iotbx.pdb.fetch.run(sys.argv[1:])
