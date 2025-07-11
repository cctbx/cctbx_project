"""
  Apply BIOMT records of PDB or equivalent records of mmCIF.
"""

from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb.biomt_reconstruction

import sys, os
import iotbx.pdb
import mmtbx.model

def run(args):
  """
  Apply BIOMT records of PDB or equivalent records of mmCIF.
  Example:
    phenix.pdb.biomt_reconstruction model.pdb
    phenix.pdb.biomt_reconstruction model.cif
  """
  if(len(args)==0 or "--help" in args or "-h" in args):
    print(run.__doc__)
    return
  file_name = args[0]
  pdb_inp = iotbx.pdb.input(file_name=file_name)
  m = mmtbx.model.manager(model_input = pdb_inp, expand_with_mtrix=False)
  m.expand_with_BIOMT_records()
  ofn = "%s_BIOMT_expanded" % (os.path.splitext(os.path.basename(file_name))[0])
  if m.input_model_format_cif():
    ofn += ".cif"
    text_to_write = m.model_as_mmcif()
  else:
    ofn += ".pdb"
    text_to_write = m.model_as_pdb()
  print("Writing result to %s file."%ofn)
  with open(ofn, 'w') as f:
    f.write(text_to_write)
if (__name__ == "__main__"):
  run(args=sys.argv[1:])
