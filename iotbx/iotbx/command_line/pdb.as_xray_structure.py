#! /usr/bin/env python

from iotbx import pdb
from iotbx.option_parser import iotbx_option_parser
from scitbx.python_utils import easy_pickle
import sys, os

def run(args):
  command_line = (iotbx_option_parser(
    usage="iotbx.pdb.as_xray_structure [options] pdb_file ...",
    description="Example: iotbx.pdb.as_xray_structure pdb1ab1.ent")
    .enable_symmetry_comprehensive()
    .option(None, "--force_symmetry",
      action="store_true",
      default=00000,
      dest="force_symmetry",
      help="symmetry on command line is stronger than symmetry found in files")
    .option("-v", "--verbose",
      action="store_true",
      default=00000,
      dest="verbose",
      help="show scatterers")
    .option(None, "--pickle",
      action="store",
      type="string",
      dest="pickle",
      help="write all data to FILE ('--pickle .' copies name of input file)",
      metavar="FILE")
  ).process(args=args)
  if (len(command_line.args) == 0):
    command_line.parser.show_help()
  all_structures = []
  for file_name in command_line.args:
    print "file_name:", file_name
    sys.stdout.flush()
    structure = pdb.as_xray_structure(
      file_name=file_name,
      crystal_symmetry=command_line.symmetry,
      force_symmetry=command_line.options.force_symmetry)
    structure.show_summary()
    if (command_line.options.verbose):
      structure.show_scatterers()
    all_structures.append(structure)
    print
  pickle_file_name = command_line.options.pickle
  if (pickle_file_name is not None and len(all_structures) > 0):
    if (pickle_file_name == "."):
      if (len(command_line.args) > 1):
        raise UserError(
          "Ambiguous name for pickle file (more than one input file).")
      pickle_file_name = os.path.basename(command_line.args[0])
    if (not pickle_file_name.lower().endswith(".pickle")):
      pickle_file_name += ".pickle"
    if (len(all_structures) == 1):
      all_structures = all_structures[0]
    print
    print "Writing all xray structures to file:", pickle_file_name
    easy_pickle.dump(pickle_file_name, all_structures)
    print

if (__name__ == "__main__"):
  run(sys.argv[1:])
