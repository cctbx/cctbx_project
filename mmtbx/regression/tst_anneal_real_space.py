
from __future__ import division
from mmtbx.regression import tst_build_alt_confs
from mmtbx.command_line import anneal_real_space
from libtbx.utils import null_out
import os.path

def exercise () :
  tst_build_alt_confs.prepare_inputs(prefix="tst_anneal_real_space")
  args = [
    "tst_anneal_real_space_start.pdb",
    "tst_anneal_real_space.mtz",
    "selection=\"chain A and resseq 2:4\"",
    "output.file_name=tst_anneal_real_space_out.pdb",
    "map_file_name=tst_anneal_real_space_maps.mtz",
    "include_starting_model=False",
  ]
  print args
  rmsd = anneal_real_space.run(args=args, out=null_out())
  assert os.path.isfile("tst_anneal_real_space_out.pdb")
  assert (rmsd > 0) and (rmsd < 1)

if (__name__ == "__main__") :
  exercise()
  print "OK"
