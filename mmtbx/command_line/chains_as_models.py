from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.chains_as_models

import sys
import iotbx.pdb
from libtbx.utils import Sorry

legend = """phenix.models_as_chains:
  Convert multi-chain PDB file into multi-model (MODEL-ENDMDL) PDB file.

How to run:
  phenix.chains_as_models model.pdb

Feedback:
  PAfonine@lbl.gov
  phenixbb@phenix-online.org"""

def run(args):
  if(len(args)!=1): raise Sorry("PDB file is expected.")
  try:
    pdb_inp = iotbx.pdb.input(file_name=args[0])
  except Exception:
    raise Sorry("PDB file is expected.")
  #
  root = iotbx.pdb.hierarchy.root()
  h = pdb_inp.construct_hierarchy()
  n_models = len(h.models())
  if(not n_models==1):
    raise Sorry("PDB file must contain one MODEL, found %d models."%n_models)
  for i, chain in enumerate(h.chains()):
    m = iotbx.pdb.hierarchy.model(id=str(i))
    chain_ = chain.detached_copy()
    chain_.id="A"
    m.append_chain(chain_)
    root.append_model(m)
  root.write_pdb_file(
    file_name        = "models_"+args[0],
    crystal_symmetry = pdb_inp.crystal_symmetry())


if (__name__ == "__main__"):
  run(args=sys.argv[1:])
