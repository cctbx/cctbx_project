import iotbx.pdb.xray_structure
from cctbx import sgtbx
from cctbx.development import random_structure
from cctbx.array_family import flex
from cStringIO import StringIO
import sys

def exercise_xray_structure(anisotropic_flag, verbose=0):
  structure = random_structure.xray_structure(
    space_group_info=sgtbx.space_group_info("P 31"),
    elements=["N","C","C","O","Si"]*2,
    volume_per_atom=500,
    min_distance=2.,
    general_positions_only=00000,
    random_u_iso=0001,
    anisotropic_flag=anisotropic_flag)
  f_abs = abs(structure.structure_factors(
    anomalous_flag=00000, d_min=2, algorithm="direct").f_calc())
  for fractional_coordinates in (00000, 0001):
    pdb_file = structure.as_pdb_file(
      remark="Title", remarks=["Any", "Thing"],
      fractional_coordinates=fractional_coordinates)
    if (0 or verbose):
      sys.stdout.write(pdb_file)
    structure_read = iotbx.pdb.as_xray_structure(
      file_iterator=StringIO(pdb_file),
      fractional_coordinates=fractional_coordinates)
    f_read = abs(f_abs.structure_factors_from_scatterers(
      xray_structure=structure_read, algorithm="direct").f_calc())
    regression = flex.linear_regression(f_abs.data(), f_read.data())
    assert regression.is_well_defined()
    if (0 or verbose):
      regression.show_summary()
    assert abs(regression.slope()-1) < 1.e-2
    assert abs(regression.y_intercept()) < 0.1

def run():
  verbose = "--Verbose" in sys.argv[1:]
  for anisotropic_flag in (00000, 0001):
    exercise_xray_structure(anisotropic_flag, verbose=verbose)
  print "OK"

if (__name__ == "__main__"):
  run()
