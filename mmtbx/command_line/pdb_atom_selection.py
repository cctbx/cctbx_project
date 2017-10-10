from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb_atom_selection

import argparse
from cctbx.array_family import flex
from libtbx.utils import plural_s
from libtbx.str_utils import show_string
import libtbx.load_env
import mmtbx.model
import iotbx.pdb
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
  pdb_inp = iotbx.pdb.input(file_name=co.file_name[0])
  model = mmtbx.model.manager(
      model_input=pdb_inp,
      process_input=True)
  atoms = model.get_atoms()
  all_bsel = flex.bool(atoms.size(), False)
  for selection_string in co.inselections:
    print(selection_string)
    isel = model.iselection(selstr=selection_string)
    all_bsel.set_selected(isel, True)
    if (not co.write_pdb_file):
      print("  %d atom%s selected" % plural_s(isel.size()))
      for atom in atoms.select(isel):
        print("    %s" % atom.format_atom_record())
  print()
  if (co.write_pdb_file):
    print("Writing file:", show_string(co.write_pdb_file))
    selected_model = model.select(all_bsel)
    if (co.cryst1_replacement_buffer_layer is not None):
      import cctbx.crystal
      crystal_symmetry = cctbx.crystal.non_crystallographic_symmetry(
        sites_cart=selected_model.get_atoms().extract_xyz(),
        buffer_layer=co.cryst1_replacement_buffer_layer)
      selected_model.set_crystal_symmetry(crystal_symmetry)
    pdb_str = selected_model.model_as_pdb()
    f = open(co.write_pdb_file, 'w')
    f.write(pdb_str)
    f.close()
    print()

if (__name__ == "__main__"):
  run(sys.argv[1:])
