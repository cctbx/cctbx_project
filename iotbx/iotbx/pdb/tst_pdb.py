import iotbx.pdb.xray_structure
from cctbx import crystal
from cctbx import sgtbx
from cctbx.development import random_structure
from cctbx.array_family import flex
import libtbx.itertbx
from cStringIO import StringIO
import sys

def exercise_format_records():
  crystal_symmetry = crystal.symmetry(
    unit_cell=(10,10,13,90,90,120),
    space_group_symbol="R 3").primitive_setting()
  assert iotbx.pdb.format_cryst1_record(crystal_symmetry=crystal_symmetry) \
    == "CRYST1    7.219    7.219    7.219  87.68  87.68  87.68 R3:R           "
  assert iotbx.pdb.format_scale_records(
    unit_cell=crystal_symmetry.unit_cell()).splitlines() \
      == ["SCALE1      0.138527 -0.005617 -0.005402        0.00000",
          "SCALE2      0.000000  0.138641 -0.005402        0.00000",
          "SCALE3      0.000000  0.000000  0.138746        0.00000"]
  assert iotbx.pdb.format_scale_records(
    fractionalization_matrix=crystal_symmetry.unit_cell()
      .fractionalization_matrix(),
    u=[-1,2,-3]).splitlines() \
      == ["SCALE1      0.138527 -0.005617 -0.005402       -1.00000",
          "SCALE2      0.000000  0.138641 -0.005402        2.00000",
          "SCALE3      0.000000  0.000000  0.138746       -3.00000"]
  assert iotbx.pdb.format_atom_record() \
      == "ATOM      0  C   DUM     1       0.000   0.000   0.000" \
        +"  1.00  0.00              "

def exercise_parser():
  for i,raw_record in enumerate("""\
LINK         O1  DDA     1                 C3  DDL     2
LINK        MN    MN   391                 OE2 GLU   217            2565
""".splitlines()):
    r = iotbx.pdb.parser.pdb_record(raw_record=raw_record)
    _1 = [r.name1,r.altLoc1,r.resName1,r.chainID1,r.resSeq1,r.iCode1,r.sym1]
    _2 = [r.name2,r.altLoc2,r.resName2,r.chainID2,r.resSeq2,r.iCode2,r.sym2]
    if (i == 0):
      assert _1 == [' O1 ', ' ', 'DDA', ' ', '   1', ' ', '      ']
      assert _2 == [' C3 ', ' ', 'DDL', ' ', '   2', ' ', '      ']
    else:
      assert _1 == ['MN  ', ' ', ' MN', ' ', ' 391', ' ', '      ']
      assert _2 == [' OE2', ' ', 'GLU', ' ', ' 217', ' ', '  2565']

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
  for res_name in (None, "res"):
    for fractional_coordinates in (00000, 0001):
      pdb_file = structure.as_pdb_file(
        remark="Title", remarks=["Any", "Thing"],
        fractional_coordinates=fractional_coordinates,
        res_name=res_name)
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
  exercise_format_records()
  exercise_parser()
  for anisotropic_flag in (00000, 0001):
    exercise_xray_structure(anisotropic_flag, verbose=verbose)
  print "OK"

if (__name__ == "__main__"):
  run()
