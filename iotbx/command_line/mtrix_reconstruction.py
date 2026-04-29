"""Apply MTRIX records of PDB or equivalent records of mmCIF"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb.mtrix_reconstruction

import sys, os
import iotbx.pdb
import iotbx.cif
import mmtbx.model

def run(args):
  """
  Apply MTRIX records of PDB or equivalent records of mmCIF.
  Example:
    phenix.pdb.mtrix_reconstruction model.pdb
    phenix.pdb.mtrix_reconstruction model.cif
  """
  if(len(args)==0 or "--help" in args or "-h" in args):
    print(run.__doc__)
    return
  file_name = args[0]
  pdb_inp = iotbx.pdb.input(file_name=file_name)
  model = mmtbx.model.manager(
      model_input=pdb_inp)
  if model.input_model_format_cif():
    out_text = model.model_as_mmcif()
    ext = ".cif"
  else:
    out_text = model.model_as_pdb()
    ext = ".pdb"
  ofn = "%s_MTRIX_expanded%s"%(
    os.path.splitext(os.path.basename(file_name))[0], ext)
  print("Writing result to %s file."%ofn)
  with open(ofn, 'w') as f:
    f.write(out_text)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
