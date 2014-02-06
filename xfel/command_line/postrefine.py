# LIBTBX_SET_DISPATCHER_NAME cxi.postrefine
from __future__ import division
from libtbx.utils import multi_out
import sys

def run(args):
  from xfel.cxi.postrefine import start
  start(args) 
  import pdb; pdb.set_trace()

if (__name__ == "__main__"):
  import sys
  run(args = sys.argv[1:])
