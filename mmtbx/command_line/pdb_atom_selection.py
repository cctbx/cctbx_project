from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb_atom_selection

from mmtbx.monomer_library import pdb_interpretation
from mmtbx.monomer_library import server
import argparse
from cctbx.array_family import flex
from libtbx.utils import plural_s
from libtbx.str_utils import show_string
import libtbx.load_env
from iotbx.pdb import write_whole_pdb_file
import sys

def run(args, command_name=libtbx.env.dispatcher_name):
  parser = argparse.ArgumentParser(
      prog=command_name,
      usage='%s pdb_file "atom_selection" [...]' % command_name)
  parser.add_argument(
      "file_name",
      nargs=1,
      help="File name of the model file")
  parser.add_argument(
      "inselections",
      help="Atom selection strings",
      nargs='+',
      )
  parser.add_argument(
      "--write-pdb-file",
      action="store",
      help="write selected atoms to new PDB file",
      default=None)
  parser.add_argument(
      "--cryst1-replacement-buffer-layer",
      action="store",
      type=float,
      help="replace CRYST1 with pseudo unit cell covering the selected"
        " atoms plus a surrounding buffer layer",
      default=None)
  co = parser.parse_args(args)
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=co.file_name[0],
    log=sys.stdout)
  print
  acp = processed_pdb_file.all_chain_proxies

  selection_cache = acp.pdb_hierarchy.atom_selection_cache()
  atoms = acp.pdb_atoms
  all_bsel = flex.bool(atoms.size(), False)
  for selection_string in co.inselections:
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
    write_whole_pdb_file(
        file_name=co.write_pdb_file,
        processed_pdb_file=processed_pdb_file,
        pdb_hierarchy=sel_hierarchy,
        crystal_symmetry=crystal_symmetry)
    print

if (__name__ == "__main__"):
  run(sys.argv[1:])
