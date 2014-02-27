
from __future__ import division
from mmtbx.regression import tst_build_alt_confs
from mmtbx.command_line import generate_disorder
from libtbx.utils import null_out
import os.path

def exercise () :
  tst_build_alt_confs.prepare_inputs(prefix="tst_generate_disorder")
  args = [
    "tst_generate_disorder_start.pdb",
    "tst_generate_disorder.mtz",
    "selection=\"chain A and resseq 2:4\"",
    "output.file_name=tst_generate_disorder_out.pdb",
    "map_file_name=tst_generate_disorder_maps.mtz",
    "include_starting_model=False",
    "random_seed=12345",
  ]
  rmsds = generate_disorder.run(args=args, out=null_out())
  assert os.path.isfile("tst_generate_disorder_out.pdb")
  # TODO rmsd ~= 1.173 on my Mac, which means it's finding the second
  # conforamtion - need to check other systems and lock this in
  print rmsds[0]
  assert (rmsds[0] > 0) and (rmsds[0] < 2)

if (__name__ == "__main__") :
  exercise()
  print "OK"
