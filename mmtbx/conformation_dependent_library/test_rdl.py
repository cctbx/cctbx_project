from __future__ import division
import sys
import time

from mmtbx.conformation_dependent_library.rdl_database import rdl_database
from mmtbx.conformation_dependent_library import rotamers as rdl

from iotbx import pdb

def get_geometry_restraints_manager(pdb_filename,
                                    #pdb_inp,
                                    #pdb_hierarchy,
                                    ):
  t0=time.time()
  from mmtbx.monomer_library import server
  from mmtbx.monomer_library import pdb_interpretation
  #lines = pdb_hierarchy.as_pdb_string(
  #  crystal_symmetry=pdb_inp.crystal_symmetry(),
  #  )
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  processed_pdb = pdb_interpretation.process(
    mon_lib_srv,
    ener_lib,
    #raw_records=lines,
    file_name=pdb_filename,
    )
  geometry_restraints_manager = processed_pdb.geometry_restraints_manager()
  print 'time',time.time()-t0
  return geometry_restraints_manager

def run(filename):
  print filename
  if 1:
    print "RDL"
    for aa in sorted(rdl_database):
      print "  %s" % aa
      for key, value in rdl_database[aa].items():
        print "    %s" % key
        for names, values in rdl_database[aa][key].items():
          print "      %s : %s" % (names, values)
    assert 0
  #
  pdb_inp = pdb.input(filename)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  geometry_restraints_manager = get_geometry_restraints_manager(filename,
                                                                #pdb_inp,
                                                                #pdb_hierarchy,
                                                                )
  pdb_hierarchy.reset_i_seq_if_necessary()
  forward, reverse = rdl.adjust_rotomer_restraints(
    pdb_hierarchy,
    geometry_restraints_manager,
    verbose=True,
    )
  rdl.adjust_rotomer_restraints(pdb_hierarchy,
                                geometry_restraints_manager,
                                i_seqs_restraints=reverse,
    )

if __name__=="__main__":
  run(sys.argv[1])
