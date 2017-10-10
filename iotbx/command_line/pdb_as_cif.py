# LIBTBX_SET_DISPATCHER_NAME phenix.pdb_as_cif

from __future__ import division
from __future__ import print_function
from builtins import str
import iotbx.pdb
import os
import sys

def run(args, out=sys.stdout):
  for file_name in args:
    try:
      assert os.path.exists(file_name)
      print("Converting %s to mmCIF format." %file_name, file=out)
      pdb_input = iotbx.pdb.pdb_input(file_name=file_name)
      hierarchy = pdb_input.construct_hierarchy(sort_atoms=False)
      pdb_atoms = hierarchy.atoms()
      pdb_atoms.set_chemical_element_simple_if_necessary()
      elements = pdb_atoms.extract_element().strip()
      if (not elements.all_ne("")) :
        n_missing = elements.count("")
        raise RuntimeError("Missing element symbol for %d atoms." % n_missing)
      basename = os.path.splitext(os.path.basename(file_name))[0]
      iotbx.cif.write_whole_cif_file(
          file_name=basename+".cif",
          pdb_hierarchy=hierarchy,
          crystal_symmetry=pdb_input.crystal_symmetry(),
          ss_annotation=pdb_input.extract_secondary_structure(),
          cif_block_name=basename,
          )
      print("  wrote %s.cif" % basename, file=out)
    # except IOError, e: # debugging variant, to see traceback
    except Exception as e:
      print("Error converting %s to mmCIF format:" %file_name, file=out)
      print(" ", str(e), file=out)
      continue

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
