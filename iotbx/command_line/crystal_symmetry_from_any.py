from iotbx import crystal_symmetry_from_any
from iotbx.pdb import format_cryst1_and_scale_records
from iotbx.cns.crystal_symmetry_utils import crystal_symmetry_as_cns_inp_defines
from iotbx import format
import sys

def run(args):
  assert len(args) == 1
  crystal_symmetry = crystal_symmetry_from_any.extract_from(args[0])
  if (crystal_symmetry is None):
    raise RuntimeError, \
      "Unknown file format or unit cell and/or space group missing from file."
  format.crystal_symmetry(crystal_symmetry)
  print
  print "\n".join(
    crystal_symmetry_as_cns_inp_defines(crystal_symmetry=crystal_symmetry))
  print
  print format_cryst1_and_scale_records(
    crystal_symmetry=crystal_symmetry,
    write_scale_records=True)
  print

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
