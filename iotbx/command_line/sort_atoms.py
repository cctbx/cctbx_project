# LIBTBX_SET_DISPATCHER_NAME iotbx.pdb.sort_atoms

from __future__ import division
from libtbx.utils import Usage
import sys
import iotbx.pdb
from iotbx.pdb import write_whole_pdb_file

master_phil_str = """
file_name = None
  .type = path
  .multiple = False
  .optional = False
  .style = hidden
"""

def show_usage():
  help_msg = """\
iotbx.pdb.sort_atoms model.pdb

Sort atoms in residues so they will be in the same order in all residues.
Also renumbers atoms (atom serial number field 7-11 columns)."""

  raise Usage(help_msg)

def run(args):
  if len(args) == 0:
    show_usage()
    return
  inp_fn = args[0]
  pdb_input = iotbx.pdb.input(
      file_name=inp_fn,
      source_info=None,
      raise_sorry_if_format_error=True)
  pdb_h = pdb_input.construct_hierarchy()
  pdb_h.sort_atoms_in_place()

  out_fn = inp_fn.strip(".pdb") + "_sorted.pdb"
  write_whole_pdb_file(
      file_name=out_fn,
      output_file=None,
      processed_pdb_file=None,
      pdb_hierarchy=pdb_h,
      crystal_symmetry=pdb_input.crystal_symmetry(),
      ss_annotation=pdb_input.extract_secondary_structure(),
      atoms_reset_serial_first_value=None,
      link_records=None,
      )


if (__name__ == "__main__") :
  run(sys.argv[1:])
