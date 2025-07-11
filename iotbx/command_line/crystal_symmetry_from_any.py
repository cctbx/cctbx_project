"""Extract crystal_symmetry from any suitable file"""
from __future__ import absolute_import, division, print_function
from iotbx import crystal_symmetry_from_any
from iotbx.pdb import format_cryst1_and_scale_records
from iotbx.cns.crystal_symmetry_utils import crystal_symmetry_as_cns_inp_defines
from iotbx import format
from libtbx.option_parser import option_parser
import libtbx.load_env
import sys

def run(args):
  if (len(args) == 0): args = ["--help"]
  command_line = (option_parser(
    usage="%s [OPTIONS] FILE..." % libtbx.env.dispatcher_name)
    .option(None, "--niggli_cell",
      action="store_true")
  ).process(args=args)
  if (len(command_line.args) == 0):
    command_line.parser.show_help()
    return
  co = command_line.options
  for arg in command_line.args:
    crystal_symmetry = crystal_symmetry_from_any.extract_from(arg)
    if (crystal_symmetry is None):
      raise RuntimeError("Unknown file format or unit cell and space group missing from file.")
    if (co.niggli_cell
          and crystal_symmetry.unit_cell() is not None
          and crystal_symmetry.space_group_info() is not None):
      crystal_symmetry = crystal_symmetry.niggli_cell()
    format.crystal_symmetry(crystal_symmetry)
    print()
    print("\n".join(
      crystal_symmetry_as_cns_inp_defines(crystal_symmetry=crystal_symmetry)))
    print()
    print(format_cryst1_and_scale_records(
      crystal_symmetry=crystal_symmetry,
      write_scale_records=True))
    print()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
