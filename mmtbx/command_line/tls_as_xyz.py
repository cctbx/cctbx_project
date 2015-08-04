from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.tls_as_xyz

import sys
from libtbx.utils import Sorry
import mmtbx.tls.tls_as_xyz

legend = """phenix.tls_as_xyz:
  Given PDB file with TLS records generate ensemble of models (multi-model PDB
  file) that are consistent with TLS parameters.

How to run:
  phenix.tls_as_xyz model.pdb

"""

def run(args, log=sys.stdout):
  print >> log, "-"*79
  print >> log, legend
  print >> log, "-"*79
  inputs = mmtbx.utils.process_command_line_args(args = args)
  file_names = inputs.pdb_file_names
  if(len(file_names) != 1): raise Sorry("PDB file must be provided.")
  mmtbx.tls.tls_as_xyz.run(pdb_file_name = file_names[0], n_models=499, log=log)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
