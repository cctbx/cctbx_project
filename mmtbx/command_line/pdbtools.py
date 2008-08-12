# LIBTBX_SET_DISPATCHER_NAME phenix.pdbtools

from mmtbx import pdbtools
import sys

if (__name__ == "__main__"):
  pdbtools.run(args=sys.argv[1:])
