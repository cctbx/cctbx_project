"""Try to convert an mmCIF formatted model file to PDB format"""

from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.cif_as_pdb
# LIBTBX_SET_DISPATCHER_NAME iotbx.cif_as_pdb
import os
import iotbx.pdb
import iotbx.pdb.mmcif
import mmtbx.model


def run(args):
  for file_name in args:
    try:
      assert os.path.exists(file_name)
      print("Converting %s to PDB format." %file_name)
      cif_input = iotbx.pdb.mmcif.cif_input(file_name=file_name)
      m = mmtbx.model.manager(model_input=cif_input)
      basename = os.path.splitext(os.path.basename(file_name))[0]
      if (not m.can_be_output_as_pdb()):
        # Convert the hierarchy to forward_compatible:
        from iotbx.pdb.forward_compatible_pdb_cif_conversion import \
           hierarchy_as_forward_compatible_pdb_string
        pdb_text =  \
              hierarchy_as_forward_compatible_pdb_string(m.get_hierarchy())
        print("***Warning: the file %s does not fit in standard PDB format." %(
         file_name) + \
         "\nConverting to forward_compatible PDB, potentially changing \nsome "+
         "chain ID and residue names")
      else:
        pdb_text = m.model_as_pdb()
      print("Writing %s" % (basename+".pdb"))
      with open(basename+".pdb", 'w') as f:
        f.write(pdb_text)
    except Exception as e:
      print("Error converting %s to PDB format:" %file_name)
      print(" ", str(e))
      continue

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
