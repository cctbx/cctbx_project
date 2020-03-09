from __future__ import division
import time

def get_geometry_restraints_manager(pdb_filename=None,
                                    raw_records=None,
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
    raw_records=raw_records,
    file_name=pdb_filename,
    )
  geometry_restraints_manager = processed_pdb.geometry_restraints_manager()
  print('time',time.time()-t0)
  return geometry_restraints_manager
