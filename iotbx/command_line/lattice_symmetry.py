"""Analyze lattice symmetry"""

from __future__ import absolute_import, division, print_function
# Comments by Phil Evans, MRC-LMB, Cambridge, U.K.

from cctbx import crystal
from cctbx.sgtbx.lattice_symmetry import metric_subgroups
from iotbx.option_parser import option_parser
import sys

def run(args):
  command_line = (option_parser(
    usage="iotbx.lattice_symmetry [options] [centring_type_symbol]",
    description="Example: iotbx.lattice_symmetry"
               +" --unit_cell=12,12,12.1,89,90,92 F")
    .enable_symmetry_comprehensive()
    .option(None, "--delta",
      action="store",
      type="float",
      default=3.,
      help="angular tolerance in degrees")
  ).process(args=args, max_nargs=1)
  # Pick up symmetry object
  input_symmetry = command_line.symmetry
  # Check that we have what we need
  if (input_symmetry.unit_cell() is None):
    print()
    print("***********************************")
    print("Please specify unit cell parameters")
    print("***********************************")
    print()
    command_line.parser.show_help()
    return
  if (len(command_line.args) > 0):
    input_symmetry = crystal.symmetry(
      unit_cell=input_symmetry.unit_cell(),
      space_group_symbol="Hall: %s 1" % command_line.args[0])
  elif (input_symmetry.space_group_info() is None):
    input_symmetry = crystal.symmetry(
      unit_cell=input_symmetry.unit_cell(),
      space_group_symbol="P 1")
  # Do it
  groups = metric_subgroups(input_symmetry, command_line.options.delta,
    enforce_max_delta_for_generated_two_folds=True)
  groups.show()

if (__name__ == "__main__"):
  run(sys.argv[1:])
