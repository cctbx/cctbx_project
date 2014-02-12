# LIBTBX_SET_DISPATCHER_NAME cxi.find_target_uc
from __future__ import division
from libtbx.utils import multi_out
import sys

#def run(args):
#  log = open("find_target_uc.log", "w")
#  out = multi_out()
#  out.register("log", log, atexit_send_to=None)
#  out.register("stdout", sys.stdout)
#
#  print >> out, "Target unit cell and space group:"
#  print >> out, "  ", work_params.target_unit_cell
#  print >> out, "  ", work_params.target_space_group
#
def run(args):
  #-------- Get the pickles from a command line specified path arg -------#
  if args < 2: raise IOError("Must give at least one path to folder of pickles") 
  from cctbx.determine_unit_cell.target_uc import target
  ucs = target(args) 
  import pdb; pdb.set_trace()

if (__name__ == "__main__"):
  import sys
  result = run(args = sys.argv[1:])
