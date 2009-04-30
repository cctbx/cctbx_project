from iotbx import crystal_symmetry_from_any
from iotbx.pdb import format_cryst1_and_scale_records
from iotbx import format
import sys

def run(args):
  assert len(args) == 1
  crystal_symmetry = crystal_symmetry_from_any.extract_from(args[0])
  if (crystal_symmetry is None):
    raise RuntimeError, \
      "Unknown file format or unit cell and/or space group missing from file."
  format.crystal_symmetry(crystal_symmetry)
  print format_cryst1_and_scale_records(
    crystal_symmetry=crystal_symmetry,
    write_scale_records=True)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
