#! /usr/bin/env python

from iotbx import pdb
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from iotbx import crystal_symmetry_from_any
from scitbx.python_utils import easy_pickle
import sys

def usage():
  return (  "usage: iotbx.pdb.as_xray_structure"
          + " [--unit_cell=1,1,1,90,90,90]"
          + " [--space_group=P212121]"
          + " [--extract_symmetry=any_file_format]"
          + " [--force_symmetry]"
          + " [--pickle=file_name]"
          + " pdb_file ...")

def run(args):
  unit_cell = None
  space_group_info = None
  force_symmetry = 00000
  pickle_file_name = None
  remaining_args = []
  for arg in args:
    if (arg.startswith("--unit_cell=")):
      params = arg.split("=", 1)[1]
      unit_cell = uctbx.unit_cell(params)
    elif (arg.startswith("--space_group=")):
      symbol = arg.split("=", 1)[1]
      space_group_info = sgtbx.space_group_info(symbol=symbol)
    elif (arg.startswith("--extract_symmetry=")):
      file_name = arg.split("=", 1)[1]
      crystal_symmetry = crystal_symmetry_from_any.extract_from(file_name)
      unit_cell = crystal_symmetry.unit_cell()
      space_group_info = crystal_symmetry.space_group_info()
    elif (arg == "--force_symmetry"):
      force_symmetry = 0001
    elif (arg.startswith("--pickle=")):
      pickle_file_name = arg.split("=", 1)[1]
    elif (arg.startswith("--")):
      print usage()
      raise RuntimeError, "Unknown option: " + arg
    else:
      remaining_args.append(arg)
  args = remaining_args
  all_structures = []
  for file_name in args:
    print "file_name:", file_name
    sys.stdout.flush()
    structure = pdb.as_xray_structure(
      file_name,
      crystal_symmetry=crystal.symmetry(
        unit_cell=unit_cell,
        space_group_info=space_group_info),
      force_symmetry=force_symmetry)
    structure.show_summary()
    all_structures.append(structure)
    print
  if (pickle_file_name is not None and len(all_structures) > 0):
    if (len(all_structures) == 1):
      all_structures = all_structures[0]
    if (not pickle_file_name.lower().endswith(".pickle")):
      pickle_file_name += ".pickle"
    print "Writing all xray structures to file:", pickle_file_name
    easy_pickle.dump(pickle_file_name, all_structures)
    print

if (__name__ == "__main__"):
  run(sys.argv[1:])
