"""Generate ensemble of models consistent with TLS parameters"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.tls_as_xyz

import sys
from libtbx.utils import Sorry
import mmtbx.tls.tls_as_xyz
import iotbx.phil

legend = """phenix.tls_as_xyz:
  Given PDB file with TLS records generate ensemble of models (multi-model PDB
  file) that are consistent with TLS parameters.

How to run:
  phenix.tls_as_xyz model.pdb

"""

master_params_str = """
n_models=499
  .type=int
"""

def master_params():
  return iotbx.phil.parse(master_params_str)

def run(args, log=sys.stdout):
  print("-"*79, file=log)
  print(legend, file=log)
  print("-"*79, file=log)
  inputs = mmtbx.utils.process_command_line_args(args = args,
    master_params = master_params())
  params = inputs.params.extract()
  file_names = inputs.pdb_file_names
  if(len(file_names) != 1): raise Sorry("PDB file must be provided.")
  output_file_name_prefix = file_names[0].replace(".pdb", "")
  mmtbx.tls.tls_as_xyz.run(
    pdb_file_name    = file_names[0],
    n_models         = params.n_models,
    log              = log,
    output_file_name_prefix = output_file_name_prefix)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])

