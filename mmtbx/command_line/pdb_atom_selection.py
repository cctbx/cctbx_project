# LIBTBX_SET_DISPATCHER_NAME phenix.pdb_atom_selection

from mmtbx.monomer_library import pdb_interpretation
from mmtbx.monomer_library import server
from iotbx.option_parser import option_parser
from cctbx.array_family import flex
from libtbx.utils import plural_s
from libtbx.str_utils import show_string
import libtbx.load_env
import sys

def run(args, command_name=libtbx.env.dispatcher_name):
  if (len(args) == 0): args = ["--help"]
  command_line = (option_parser(
    usage='%s pdb_file "atom_selection" [...]' % command_name)
    .option(None, "--write_pdb_file",
      action="store",
      type="string",
      default=None,
      help="write selected atoms to new PDB file",
      metavar="FILE")
    .option(None, "--cryst1_replacement_buffer_layer",
      action="store",
      type="float",
      default=None,
      help="replace CRYST1 with pseudo unit cell covering the selected"
        " atoms plus a surrounding buffer layer",
      metavar="WIDTH")
  ).process(args=args, min_nargs=2)
  co = command_line.options
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=command_line.args[0],
    log=sys.stdout)
  print
  acp = processed_pdb_file.all_chain_proxies
  selection_cache = acp.pdb_hierarchy.atom_selection_cache()
  atoms = acp.pdb_atoms
  all_bsel = flex.bool(atoms.size(), False)
  for selection_string in command_line.args[1:]:
    print selection_string
    isel = acp.iselection(string=selection_string, cache=selection_cache)
    all_bsel.set_selected(isel, True)
    if (not co.write_pdb_file):
      print "  %d atom%s selected" % plural_s(isel.size())
      for atom in atoms.select(isel):
        print "    %s" % atom.format_atom_record()
  print
  if (co.write_pdb_file):
    print "Writing file:", show_string(co.write_pdb_file)
    sel_hierarchy = acp.pdb_hierarchy.select(all_bsel)
    if (co.cryst1_replacement_buffer_layer is None):
      crystal_symmetry = acp.special_position_settings
    else:
      import cctbx.crystal
      crystal_symmetry = cctbx.crystal.non_crystallographic_symmetry(
        sites_cart=sel_hierarchy.atoms().extract_xyz(),
        buffer_layer=co.cryst1_replacement_buffer_layer)
    sel_hierarchy.write_pdb_file(
      file_name=co.write_pdb_file,
      crystal_symmetry=crystal_symmetry,
      append_end=True)
    print

if (__name__ == "__main__"):
  run(sys.argv[1:])
