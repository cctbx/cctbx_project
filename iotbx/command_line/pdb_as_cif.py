from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb_as_cif
import os
import iotbx.cif.model
import iotbx.pdb

def run(args):
  for file_name in args:
    try:
      assert os.path.exists(file_name)
      print "Converting %s to mmCIF format." %file_name
      pdb_input = iotbx.pdb.pdb_input(file_name=file_name)
      hierarchy = pdb_input.construct_hierarchy()
      cif_object = iotbx.cif.model.cif()
      basename = os.path.splitext(os.path.basename(file_name))[0]
      cif_object[basename] = hierarchy.as_cif_block(
        crystal_symmetry=pdb_input.crystal_symmetry())
      f = open(basename+".cif", "wb")
      print >> f, cif_object
      f.close()
    except Exception, e:
      print "Error converting %s to mmCIF format:" %file_name
      print " ", str(e)
      continue

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
