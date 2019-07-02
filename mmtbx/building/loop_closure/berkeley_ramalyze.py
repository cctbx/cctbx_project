from __future__ import absolute_import, division, print_function
import iotbx.pdb
from mmtbx.building.loop_closure.utils import list_rama_outliers_h
from mmtbx.rotamer import ramachandran_eval

import sys

def run(args, log=sys.stdout):
  pdb_h = iotbx.pdb.input(source_info=None, file_name=args[0]).\
      construct_hierarchy()
  r = ramachandran_eval.RamachandranEval()
  outp = list_rama_outliers_h(pdb_h, r.rama_eval)
  print(outp)
  print("END")

if (__name__ == "__main__"):
  run(sys.argv[1:])
