from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb.mtrix_reconstruction

import sys, os
import iotbx.pdb

def run(args):
  """
  Apply MTRIX records of PDB or equivalent records of mmCIF.
  Example:
    phenix.pdb.mtrix_reconstruction model.pdb
    phenix.pdb.mtrix_reconstruction model.cif
  """
  if(len(args)==0 or "--help" in args or "-h" in args):
    print mtrix_reconstruction.__doc__
  file_name = args[0]
  pdb_inp = iotbx.pdb.input(file_name=file_name)
  h = pdb_inp.construct_hierarchy_MTRIX_expanded()
  ofn = "%s_MTRIX_expanded.pdb"%os.path.splitext(os.path.basename(file_name))[0]
  print "Writing result to %s file."%ofn
  h.write_pdb_file(file_name=ofn)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
