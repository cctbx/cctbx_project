# LIBTBX_SET_DISPATCHER_NAME phenix.secondary_structure_restraints

if __name__ == "__main__" :
  from mmtbx import secondary_structure
  import sys
  secondary_structure.run(sys.argv[1:])
