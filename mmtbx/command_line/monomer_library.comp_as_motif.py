from mmtbx.monomer_library import server
import sys

def run(args):
  list_cif = server.mon_lib_list_cif()
  srv = server.server(list_cif=list_cif)
  for comp_id in args:
    comp_comp_id = srv.get_comp_comp_id_direct(comp_id=comp_id)
    motif = comp_comp_id.as_geometry_restraints_motif()
    motif.show()

if (__name__ == "__main__"):
  run(sys.argv[1:])
