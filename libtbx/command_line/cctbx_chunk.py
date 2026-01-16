# LIBTBX_SET_DISPATCHER_NAME libtbx.chunk
from __future__ import absolute_import, division, print_function
from libtbx.queuing_system_utils import chunk_manager
from libtbx.utils import Usage
import sys

def usage():
  raise Usage("libtbx.chunk [n] [i]")

def run(args):
  if (len(args) == 0): usage()
  for arg in args:
    if (arg not in ["n", "i"]): usage()
  cm = chunk_manager(n=1, i=0).queuing_system_overrides_chunk()
  for arg in args:
    if   (arg == "n"): print(cm.n, end=' ')
    elif (arg == "i"): print(cm.i, end=' ')
  print()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
