from iotbx import pdb
import iotbx.pdb.remark_290_interpretation
import iotbx.pdb.atom
import iotbx.pdb.interpretation
import iotbx.pdb.xray_structure
from cctbx import crystal
from cctbx import sgtbx
from cctbx.development import random_structure
from cctbx.array_family import flex
import libtbx.itertbx
from cStringIO import StringIO
import sys, os

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
      assert _1 == [' O1 ', ' ', 'DDA', ' ', 1, ' ', '      ']
      assert _2 == [' C3 ', ' ', 'DDL', ' ', 2, ' ', '      ']
    else:
      assert _1 == ['MN  ', ' ', ' MN', ' ', 391, ' ', '      ']
      assert _2 == [' OE2', ' ', 'GLU', ' ', 217, ' ', '  2565']

def exercise_remark_290_interpretation():
  symmetry_operators=pdb.remark_290_interpretation.extract_symmetry_operators(
    remark_290_records=pdb.remark_290_interpretation.example.splitlines())
  assert symmetry_operators is not None
  assert len(symmetry_operators) == 4
  assert symmetry_operators[0] == sgtbx.rt_mx("X,Y,Z")
  assert symmetry_operators[1] == sgtbx.rt_mx("1/2-X,-Y,1/2+Z")
  assert symmetry_operators[2] == sgtbx.rt_mx("-X,1/2+Y,1/2-Z")
  assert symmetry_operators[3] == sgtbx.rt_mx("1/2+X,1/2-Y,-Z")
  for link_sym,expected_sym_op in [("1555", "x,y,z"),
                                   ("1381", "x-2,y+3,z-4"),
                                   ("3729", "-x+2,1/2+y-3,1/2-z+4"),
                                   (" 3_729 ", "-x+2,1/2+y-3,1/2-z+4"),
                                   (" 3 729 ", "-x+2,1/2+y-3,1/2-z+4"),
                                   ("_3729", None),
                                   ("37_29", None)]:
    sym_op = pdb.remark_290_interpretation.get_link_symmetry_operator(
      symmetry_operators=symmetry_operators,
      link_sym=link_sym)
    if (sym_op is None):
      assert expected_sym_op is None
    else:
      assert sym_op == sgtbx.rt_mx(expected_sym_op)

def exercise_atom():
  for lbls in [pdb.atom.labels(),
               pdb.atom.labels_from_string(" C  ,,,,,,,"),
               pdb.atom.labels_from_string(" C  ,A,GLY,C,1,I,S,4")]:
    assert str(lbls) == str(pdb.atom.labels_from_string(str(lbls)))
  atom_etc_records = """\
ATOM     10  N  ACYS "   6       4.535   4.190  -0.757  1.00  0.05      1ETN N1+
SIGATM   10  N  ACYS "   6       0.010   0.210   0.310  0.02  0.03      1ETN N1+
ANISOU   10  N  ACYS "   6      441    580    402     48    -72    -88  1ETN N1+
SIGUIJ   10  N  ACYS "   6        3      2      8      3      8      6  1ETN N1+
""".splitlines()
  attr = pdb.atom.attributes()
  attr.set_from_ATOM_record(pdb.parser.pdb_record(atom_etc_records[0]))
  attr.set_from_SIGATM_record(pdb.parser.pdb_record(atom_etc_records[1]))
  attr.set_from_ANISOU_record(pdb.parser.pdb_record(atom_etc_records[2]))
  attr.set_from_SIGUIJ_record(pdb.parser.pdb_record(atom_etc_records[3]))
  assert str(attr) == " N  ,A,CYS,\",6, ,1ETN,None"
  s = StringIO()
  attr.show(f=s, prefix="  ")
  assert s.getvalue() == """\
  name:        " N  "
  altLoc:      "A"
  resName:     "CYS"
  chainID:     "\\""
  resSeq:      6
  iCode:       " "
  segID:       "1ETN"
  element:     " N"
  charge:      "1+"
  coordinates: (4.535, 4.19, -0.757)
  sigCoor:     (0.01, 0.21, 0.31)
  occupancy:   1
  sigOcc:      0.02
  tempFactor:  0.05
  sigTemp:     0.03
  Ucart:       (0.0441, 0.058, 0.0402, 0.0048, -0.0072, -0.0088)
  sigUcart:    (0.0003, 0.0002, 0.0008, 0.0003, 0.0008, 0.0006)
"""

def exercise_altLoc_grouping():
  altLoc_groups = pdb.interpretation.altLoc_grouping()
  altLoc_groups.add_group(["A", " "])
  assert altLoc_groups.group_list == [[" ", "A"]]
  altLoc_groups.add_group(["B", "C"])
  assert altLoc_groups.group_list == [[" ", "A"], ["B", "C"]]
  altLoc_groups.add_group(["D", "C"])
  assert altLoc_groups.group_list == [[" ", "A"], ["B", "C", "D"]]
  replacements = altLoc_groups.get_false_blank_altLoc_replacements()
  assert len(replacements) == 1
  assert replacements[0] == "X"
  altLoc_groups.add_group(["B", "A", "E"])
  assert altLoc_groups.group_list == [[" ", "A", "B", "C", "D", "E"]]
  altLoc_groups.add_group(["1", " ", "3"])
  assert altLoc_groups.group_list == [[" ", "A", "B", "C", "D", "E"],
                                      [" ", "1", "3"]]
  replacements = altLoc_groups.get_false_blank_altLoc_replacements()
  assert len(replacements) == 2
  assert replacements[0] == "X"
  assert replacements[1] == "0"
  altLoc_groups.add_group(["Q", "R"])
  assert len(altLoc_groups.group_list) == 3
  altLoc_groups.add_group(["Q", "B", "1"])
  assert altLoc_groups.group_list \
    == [[" ", "1", "3", "A", "B", "C", "D", "E", "Q", "R"]]
  replacements = altLoc_groups.get_false_blank_altLoc_replacements()
  assert len(replacements) == 1
  assert replacements[0] == "X"

def exercise_interpretation(verbose=0, quick=0001):
  pdb_dir = os.path.join(os.environ["LIBTBX_DIST_ROOT"],
    "regression", "pdb")
  if (not os.path.isdir(pdb_dir)): return
  for file_name in os.listdir(pdb_dir):
    if (    not file_name.endswith(".pdb")
        and not file_name.endswith(".ent")): continue
    if (0 or verbose):
      print file_name
    stage_1 = pdb.interpretation.stage_1(
      file_name=os.path.join(pdb_dir, file_name))
    sel_cache = stage_1.selection_cache()
    isel_name = sel_cache.get_name(" CA ")
    isel_altLoc = sel_cache.get_altLoc("A")
    isel_resName = sel_cache.get_resName("?L?")
    isel_chainID = sel_cache.get_chainID("A")
    isel_resSeq = sel_cache.get_resSeq(1)
    isel_iCode = sel_cache.get_iCode("A")
    isel_segID = sel_cache.get_segID("    ")
    isel_MODELserial = sel_cache.get_MODELserial(1)
    isel_element = sel_cache.get_element(" C")
    isel_charge = sel_cache.get_charge("2+")
    if (0 or verbose):
      print "  n_seq:", len(stage_1.atom_attributes_list)
      print "  isel_name:", [isel.size() for isel in isel_name]
      print "  isel_altLoc:", [isel.size() for isel in isel_altLoc]
      print "  isel_resName:", [isel.size() for isel in isel_resName]
      print "  isel_chainID:", [isel.size() for isel in isel_chainID]
      print "  isel_resSeq:", [isel.size() for isel in isel_resSeq]
      print "  isel_iCode:", [isel.size() for isel in isel_iCode]
      print "  isel_segID:", [isel.size() for isel in isel_segID]
      print "  isel_MODELserial:", [isel.size() for isel in isel_MODELserial]
      print "  isel_element:", [isel.size() for isel in isel_element]
      print "  isel_charge:", [isel.size() for isel in isel_charge]
    sel_name = sel_cache.sel_name(" CA ")
    sel_altLoc = sel_cache.sel_altLoc("A")
    sel_resName = sel_cache.sel_resName("?L?")
    sel_chainID = sel_cache.sel_chainID("A")
    sel_resSeq = sel_cache.sel_resSeq(1)
    sel_iCode = sel_cache.sel_iCode("A")
    sel_segID = sel_cache.sel_segID("    ")
    sel_MODELserial = sel_cache.sel_MODELserial(1)
    sel_element = sel_cache.sel_element(" C")
    sel_charge = sel_cache.sel_charge("2+")
    if (0 or verbose):
      print "  sel_name:", sel_name.count(0001)
      print "  sel_altLoc:", sel_altLoc.count(0001)
      print "  sel_resName:", sel_resName.count(0001)
      print "  sel_chainID:", sel_chainID.count(0001)
      print "  sel_resSeq:", sel_resSeq.count(0001)
      print "  sel_iCode:", sel_iCode.count(0001)
      print "  sel_segID:", sel_segID.count(0001)
      print "  sel_MODELserial:", sel_MODELserial.count(0001)
      print "  sel_element:", sel_element.count(0001)
      print "  sel_charge:", sel_charge.count(0001)
    if (quick): break

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
  verbose = "--verbose" in sys.argv[1:]
  exercise_format_records()
  exercise_remark_290_interpretation()
  exercise_parser()
  exercise_atom()
  exercise_altLoc_grouping()
  exercise_interpretation(verbose=verbose)
  for anisotropic_flag in (00000, 0001):
    exercise_xray_structure(anisotropic_flag, verbose=verbose)
  print "OK"

if (__name__ == "__main__"):
  run()
