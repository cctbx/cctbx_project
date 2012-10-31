# LIBTBX_SET_DISPATCHER_NAME mmtbx.find_torsion_angle_ncs_groups
from __future__ import division
import sys
import time
from mmtbx.torsion_restraints import torsion_ncs
from iotbx import pdb
from libtbx.utils import Sorry

def run(args):
  if len(args) != 1:
    raise Sorry("mmtbx.find_torsion_angle_ncs_groups requires "+
                "one PDB files as input")
  file_name = args[0]
  pdb_io = pdb.input(file_name)
  pdb_hierarchy = pdb_io.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  ncs_groups = torsion_ncs.determine_ncs_groups(
                 pdb_hierarchy=pdb_hierarchy)
  if len(ncs_groups) == 0:
    print "No NCS groups found"
  for i, group in enumerate(ncs_groups):
    print "Group %d:" % (i+1)
    for chain in group:
      print "  "+chain
    print

if (__name__ == "__main__"):
  t0 = time.time()
  run(args=sys.argv[1:])
  print "Time: %10.3f"%(time.time()-t0)
