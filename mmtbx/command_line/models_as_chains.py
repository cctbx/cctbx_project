from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.models_as_chains

import sys
import iotbx.pdb
from libtbx.utils import Sorry
import string

legend = """phenix.models_as_chains:
  Convert multi-model PDB file (MODEL-ENDMDL) into multi-chain PDB file.

How to run:
  phenix.models_as_chains model.pdb

Feedback:
  PAfonine@lbl.gov
  phenixbb@phenix-online.org"""

def run(args):
  if(len(args)!=1): raise Sorry("PDB file is expected.")
  try:
    pdb_inp = iotbx.pdb.input(file_name=args[0])
  except Exception:
    raise Sorry("PDB file is expected.")
  h = pdb_inp.construct_hierarchy()
  r = iotbx.pdb.hierarchy.root()
  m = iotbx.pdb.hierarchy.model()
  idl = [i for i in string.ascii_lowercase]
  idu = [i for i in string.ascii_uppercase]
  taken = []
  c1 = None
  c2 = None
  n_atoms = []
  for m_ in h.models():
    for c_ in m_.chains():
      n_at = len(c_.atoms())
      if(not n_at in n_atoms): n_atoms.append(n_at)
      c_ = c_.detached_copy()
      found = False
      for idu_ in idu:
        for idl_ in idl:
          id_ = idu_+idl_
          if(not id_ in taken):
            taken.append(id_)
            found = id_
            break
        if(found): break
      c_.id = found
      m.append_chain(c_)
  r.append_model(m)
  r.write_pdb_file(
    file_name        = "chains_"+args[0],
    crystal_symmetry = pdb_inp.crystal_symmetry())

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
