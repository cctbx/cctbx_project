# LIBTBX_SET_DISPATCHER_NAME phenix.pdb_as_cif

from __future__ import division
import iotbx.cif.model
import iotbx.pdb
import os
import sys

def run(args, out=sys.stdout):
  for file_name in args:
    try:
      assert os.path.exists(file_name)
      print >> out, "Converting %s to mmCIF format." %file_name
      pdb_input = iotbx.pdb.pdb_input(file_name=file_name)
      hierarchy = pdb_input.construct_hierarchy()
      pdb_atoms = hierarchy.atoms()
      pdb_atoms.set_chemical_element_simple_if_necessary()
      elements = pdb_atoms.extract_element().strip()
      if (not elements.all_ne("")) :
        n_missing = elements.count("")
        raise RuntimeError("Missing element symbol for %d atoms." % n_missing)
      cif_object = iotbx.cif.model.cif()
      basename = os.path.splitext(os.path.basename(file_name))[0]
      cif_object[basename] = hierarchy.as_cif_block(
        crystal_symmetry=pdb_input.crystal_symmetry())
      f = open(basename+".cif", "wb")
      print >> f, cif_object
      print >> out, "  wrote %s.cif" % basename
      f.close()
    except Exception, e:
      print >> out, "Error converting %s to mmCIF format:" %file_name
      print >> out, " ", str(e)
      continue

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
