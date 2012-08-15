from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.tls

from mmtbx.tls import command_line
import sys

if(__name__ == "__main__"):
  command_line.run(sys.argv[1:])
