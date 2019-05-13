
from __future__ import division
from __future__ import print_function
from mmtbx.command_line import partial_omit_map
from mmtbx.regression import tst_build_alt_confs
from libtbx.utils import null_out

def exercise():
  tst_build_alt_confs.prepare_inputs("tst_partial_omit_map")
  args = [
    "tst_partial_omit_map_start.pdb",
    "tst_partial_omit_map.mtz",
    "remove_waters=False",
    "occ=0.6",
    "nproc=1",
    "ccp4_map=partial_omit_map.ccp4",
  ]
  partial_omit_map.run(args=args, out=null_out())

if (__name__ == "__main__"):
  exercise()
  print("OK")
