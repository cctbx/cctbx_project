from __future__ import absolute_import, division, print_function

import libtbx.load_env
import os.path
from libtbx import easy_run

def run():
  """
  Just make sure it runs.
  """
  pdb_fname = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/pdbs/p9.pdb", test=os.path.isfile)
  seq_filename = "tst_model_vs_sequence_seq.fa"
  assert not easy_run.call("iotbx.pdb_as_fasta %s output_file=%s" % (pdb_fname, seq_filename))
  cmd = "mmtbx.model_vs_sequence %s %s" % (pdb_fname, seq_filename)
  print(cmd)
  assert not easy_run.call(cmd)

if (__name__ == "__main__"):
  run()
