"""Try to convert an mmCIF formatted model file to PDB format"""

from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.cif_as_pdb
# LIBTBX_SET_DISPATCHER_NAME iotbx.cif_as_pdb
import os
import iotbx.pdb
import iotbx.pdb.mmcif
import mmtbx.model


def check_args(args):
  #  Look for "force_pdb_format=True" in args and set force_pdb_format if so
  force_pdb_format = False
  for arg in args:
    if arg.lower() == "force_pdb_format=true":
       force_pdb_format = True
       args.remove(arg)
  return args, force_pdb_format

def run(args):
  args, force_pdb_format = check_args(args)
  for file_name in args:
    try:
      assert os.path.exists(file_name)
      print("Converting %s to PDB format." %file_name)
      cif_input = iotbx.pdb.mmcif.cif_input(file_name=file_name)
      m = mmtbx.model.manager(model_input=cif_input)
      basename = os.path.splitext(os.path.basename(file_name))[0]
      if (not m.can_be_output_as_pdb()):
        if not force_pdb_format:
           from libtbx.utils import Sorry
           raise Sorry(
             "The file %s cannot be converted to standard PDB format." %(
             file_name) + "\n\nUse 'force_pdb_format=True' to convert using"+
             " forward-compatible PDB \nwith 2-character chain IDs and " +
             "3-character residue names")
        else: # go ahead
          # Convert the hierarchy to forward_compatible:
          from iotbx.pdb.forward_compatible_pdb_cif_conversion import \
             hierarchy_as_forward_compatible_pdb_string
          pdb_text =  \
              hierarchy_as_forward_compatible_pdb_string(m.get_hierarchy())
          print(
           "\n***Warning: the file %s does not fit in standard PDB format." %(
           file_name) + \
           "\nConverting to forward_compatible PDB %s, with" %(
             basename+".pdb") +\
            "\n2-character chain ID and 3-character residue names.\n")
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
