import iotbx.cns.xray_structure
from iotbx.cns import sdb_reader
from cctbx import crystal
from cctbx import sgtbx
from cctbx.development import random_structure
from cctbx.array_family import flex
from cStringIO import StringIO
import sys

def exercise_sdb(verbose=0):
  structure = random_structure.xray_structure(
    space_group_info=sgtbx.space_group_info("P 31"),
    elements=["N","C","C","O"]*2,
    volume_per_atom=500,
    min_distance=2.,
    general_positions_only=00000,
    random_u_iso=0001)
  f_abs = abs(structure.structure_factors(
    anomalous_flag=00000, d_min=2, algorithm="direct").f_calc())
  sdb_out = structure.as_cns_sdb_file(
    file="foo.sdb",
    description="random_structure",
    comment=["any", "thing"],
    group="best")
  if (0 or verbose):
    sys.stdout.write(sdb_out)
  sdb_files = sdb_reader.multi_sdb_parser(StringIO(sdb_out))
  assert len(sdb_files) == 1
  structure_read = sdb_files[0].as_xray_structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=structure.unit_cell(),
      space_group_info=None))
  f_read = abs(f_abs.structure_factors_from_scatterers(
    xray_structure=structure_read, algorithm="direct").f_calc())
  regression = flex.linear_regression(f_abs.data(), f_read.data())
  assert regression.is_well_defined()
  if (0 or verbose):
    regression.show_summary()
  assert abs(regression.slope()-1) < 1.e-4
  assert abs(regression.y_intercept()) < 1.e-3

def run():
  verbose = "--Verbose" in sys.argv[1:]
  exercise_sdb(verbose)
  print "OK"

if (__name__ == "__main__"):
  run()
