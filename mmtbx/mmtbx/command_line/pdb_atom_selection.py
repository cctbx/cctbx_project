# LIBTBX_SET_DISPATCHER_NAME phenix.pdb_atom_selection

from mmtbx.monomer_library import pdb_interpretation
from mmtbx.monomer_library import server
from libtbx.utils import Usage, plural_s
import libtbx.load_env
import sys, os

def run(args, command_name=libtbx.env.dispatcher_name):
  if (len(args) < 2 or not os.path.isfile(args[0])):
    raise Usage(
      '%s pdb_file "atom_selection" [...]' % libtbx.env.dispatcher_name)
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=args[0],
    log=sys.stdout)
  print
  acp = processed_pdb_file.all_chain_proxies
  selection_cache = acp.stage_1.selection_cache()
  aal = acp.stage_1.atom_attributes_list
  for selection_string in args[1:]:
    print selection_string
    isel = acp.iselection(string=selection_string, cache=selection_cache)
    print "  %d atom%s selected" % plural_s(isel.size())
    for i_seq in isel:
      print "    %s" % aal[i_seq].pdb_format()
    print

if (__name__ == "__main__"):
  run(sys.argv[1:])
