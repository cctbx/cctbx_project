from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb.biomt_reconstruction

import sys, os
import iotbx.pdb
import iotbx.cif

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
  input_file_is_cif = isinstance(pdb_inp, iotbx.pdb.mmcif.cif_input)
  h = pdb_inp.construct_hierarchy_BIOMT_expanded()
  ss_annot = pdb_inp.construct_ss_annotation_expanded(exp_type='biomt')
  if(input_file_is_cif):
    ext=".cif"
    wff = iotbx.cif.write_whole_cif_file
  else:
    ext=".pdb"
    wff = iotbx.pdb.write_whole_pdb_file
  ofn = "%s_BIOMT_expanded%s"%(
    os.path.splitext(os.path.basename(file_name))[0], ext)
  print("Writing result to %s file."%ofn)
  wff(
    file_name        = ofn,
    pdb_hierarchy    = h,
    crystal_symmetry = pdb_inp.crystal_symmetry(),
    ss_annotation    = ss_annot)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
