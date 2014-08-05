# LIBTBX_SET_DISPATCHER_NAME mmtbx.afitt

from __future__ import division
import sys
from mmtbx.geometry_restraints import afitt

if __name__ == "__main__":
  print sys.argv
  if len(sys.argv[1:])<3:
    print '''
  Usage:
    mmtbx.afitt model_file_name restraints_file_name ligand_code
    '''
  else:
    afitt.run2()
