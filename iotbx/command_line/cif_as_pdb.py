from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.cif_as_pdb
# LIBTBX_SET_DISPATCHER_NAME iotbx.cif_as_pdb
from builtins import str
import os
import iotbx.pdb
import iotbx.pdb.mmcif

def run(args):
  for file_name in args:
    try:
      assert os.path.exists(file_name)
      print("Converting %s to PDB format." %file_name)
      cif_input = iotbx.pdb.mmcif.cif_input(file_name=file_name)
      hierarchy = cif_input.construct_hierarchy()
      basename = os.path.splitext(os.path.basename(file_name))[0]
      iotbx.pdb.write_whole_pdb_file(
          file_name=basename+".pdb",
          output_file=None,
          processed_pdb_file=None,
          pdb_hierarchy=hierarchy,
          crystal_symmetry=cif_input.crystal_symmetry(),
          ss_annotation=cif_input.extract_secondary_structure(),
          append_end=True,
          atoms_reset_serial_first_value=None,
          link_records=None)
    except Exception as e:
      print("Error converting %s to PDB format:" %file_name)
      print(" ", str(e))
      continue

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
