from iotbx import pdb
import iotbx.pdb.remark_290_interpretation
import iotbx.pdb.atom
import iotbx.pdb.interpretation
import iotbx.pdb.xray_structure
from cctbx import crystal
from cctbx import sgtbx
from cctbx.development import random_structure
from cctbx.array_family import flex
from cStringIO import StringIO
import sys, os

def exercise_format_records():
  crystal_symmetry = crystal.symmetry(
    unit_cell=(10,10,13,90,90,120),
    space_group_symbol="R 3").primitive_setting()
  assert iotbx.pdb.format_cryst1_record(crystal_symmetry=crystal_symmetry) \
    == "CRYST1    7.219    7.219    7.219  87.68  87.68  87.68 R3:R"
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
      == "ATOM      0  C   DUM     1       0.000   0.000   0.000  1.00  0.00"
  assert iotbx.pdb.format_ter_record() \
      == "TER       0      DUM     1"

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
  record name: ATOM
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

def exercise_interpretation(verbose=0, quick=True):
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
      print "  sel_name:", sel_name.count(True)
      print "  sel_altLoc:", sel_altLoc.count(True)
      print "  sel_resName:", sel_resName.count(True)
      print "  sel_chainID:", sel_chainID.count(True)
      print "  sel_resSeq:", sel_resSeq.count(True)
      print "  sel_iCode:", sel_iCode.count(True)
      print "  sel_segID:", sel_segID.count(True)
      print "  sel_MODELserial:", sel_MODELserial.count(True)
      print "  sel_element:", sel_element.count(True)
      print "  sel_charge:", sel_charge.count(True)
    if (quick): break

def exercise_scalen():
  pdb_file = """\
REMARK   4 napt COMPLIES WITH FORMAT V. 2.0
CRYST1    8.098    5.953    8.652  90.00 124.40  90.00 P 21/c
SCALE1       0.14966   0.00000   0.00000        0.00000
SCALE2       0.00000   0.16798   0.00000        0.00000
SCALE3       0.07914   0.00000   0.11558        0.00000
HETATM    1  C1  UNK     0       0.550   0.110   2.464  1.00  0.00           C
HETATM    2  X2  UNK     0       0.755   0.975   1.412  1.00  0.00           C
HETATM    3 8H   UNK     0       0.947   2.326  -0.853  1.00  0.00           H
HETATM    4  C3  UNK     0       0.321   0.626   0.102  1.00  0.00           C
HETATM    5 7H   UNK     0       1.194   1.819   1.528  1.00  0.00
HETATM    6  C_4 UNK     0       0.512   1.499  -1.006  1.00  0.00
HETATM    7 6H   UNK     0       0.830   0.351   3.372  1.00  0.00
HETATM    8  C5  UNK     0      -0.088  -1.132   2.263  1.00  0.00
HETATM    9 9H   UNK     0      -0.223  -1.757   3.019  1.00  0.00
TER      10      UNK     0
END
""".splitlines()
  stage_1 = pdb.interpretation.stage_1(raw_records=pdb_file)
  xray_structure = stage_1.extract_xray_structure(
    infer_scattering_types_from_names=True)
  s = StringIO()
  xray_structure.show_summary(f=s).show_scatterers(f=s)
  assert s.getvalue() == """\
Number of scatterers: 9
At special positions: 0
Unit cell: (8.098, 5.953, 8.652, 90, 124.4, 90)
Space group: P 1 21/c 1 (No. 14)
Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso
0    C      4 ( 0.0823  0.0185  0.3283) 1.00 0.0000
1    C      4 ( 0.1130  0.1638  0.2229) 1.00 0.0000
2    H      4 ( 0.1417  0.3907 -0.0236) 1.00 0.0000
3    C      4 ( 0.0480  0.1052  0.0372) 1.00 0.0000
4    H      4 ( 0.1787  0.3056  0.2711) 1.00 0.0000
5    C      4 ( 0.0766  0.2518 -0.0758) 1.00 0.0000
6    H      4 ( 0.1242  0.0590  0.4554) 1.00 0.0000
7    C      4 (-0.0132 -0.1902  0.2546) 1.00 0.0000
8    H      4 (-0.0334 -0.2951  0.3313) 1.00 0.0000
"""

def exercise_selection():
  pdb_file = """\
CRYST1   50.800   50.800  155.300  90.00  90.00  90.00 P 43 21 2     8
MODEL        1
ATOM      4  N   SER     1       8.753  29.755  61.685  1.00 49.13
ATOM      5  CA  SER     1       9.242  30.200  62.974  1.00 46.62
ATOM      6  C   SER     1      10.453  29.500  63.579  1.00 41.99
ATOM      7  O   SER     1      10.593  29.607  64.814  1.00 43.24
ATOM      8  CB  SER     1       8.052  30.189  63.974  1.00 53.00
ATOM      9  OG  SER     1       7.294  31.409  63.930  1.00 57.79
ATOM     10  N   ARG     2      11.360  28.819  62.827  1.00 36.48
ATOM     11  CA  ARG     2      12.548  28.316  63.532  1.00 30.20
ATOM     12  C   ARG     2      13.502  29.501  63.500  1.00 25.54
ATOM     13  O   ARG     2      13.730  30.037  62.407  1.00 23.86
ATOM     14  CB  ARG     2      13.241  27.119  62.861  1.00 27.44
ATOM     15  CG  ARG     2      12.412  25.849  62.964  1.00 23.66
ATOM     16  CD  ARG     2      13.267  24.651  63.266  1.00 23.98
ATOM     17  NE  ARG     2      13.948  24.115  62.135  1.00 22.71
ATOM     18  CZ  ARG     2      15.114  23.487  62.201  1.00 21.38
ATOM     19  NH1 ARG     2      15.845  23.331  63.301  1.00 19.34
ATOM     20  NH2 ARG     2      15.575  23.030  61.051  1.00 26.66
ATOM     21  N   PRO     3J     13.947  29.997  64.680  1.00 22.94
ATOM     22  CA  PRO     3J     14.902  31.100  64.827  1.00 20.19
ATOM     23  C   PRO     3J     16.195  30.718  64.086  1.00 18.44
ATOM     24  O   PRO     3J     16.545  29.521  64.086  1.00 19.76
ATOM     25  CB  PRO     3J     15.133  31.218  66.313  1.00 19.17
ATOM     26  CG  PRO     3J     14.065  30.364  66.951  1.00 15.12
ATOM     27  CD  PRO     3J     13.816  29.289  65.966  1.00 19.56
ATOM     28  N  AILE     4      16.953  31.648  63.512  1.00 15.29
ATOM     29  CA AILE     4      18.243  31.372  62.859  1.00 14.32
ATOM     30  C  AILE     4      19.233  32.112  63.743  1.00 13.54
ATOM     31  O  AILE     4      19.105  33.315  64.009  1.00 11.84
ATOM     32  CB AILE     4      18.298  31.951  61.406  1.00 13.62
ATOM     33  CG1AILE     4      17.157  31.300  60.620  1.00 18.39
ATOM     34  CG2AILE     4      19.661  31.747  60.743  1.00 13.64
ATOM     35  CD1AILE     4      16.879  32.102  59.355  1.00 16.69
ATOM     28  N  BILE     4      16.953  31.648  63.512  1.00 15.29
ATOM     29  CA BILE     4      18.243  31.372  62.859  1.00 14.32
ATOM     30  C  BILE     4      19.233  32.112  63.743  1.00 13.54
ATOM     31  O  BILE     4      19.105  33.315  64.009  1.00 11.84
ATOM     32  CB BILE     4      18.298  31.951  61.406  1.00 13.62
ATOM     33  CG1BILE     4      17.157  31.300  60.620  1.00 18.39
ATOM     34  CG2BILE     4      19.661  31.747  60.743  1.00 13.64
ATOM     35  CD1BILE     4      16.879  32.102  59.355  1.00 16.69
TER      36      ILE     4
ENDMDL
MODEL        2
HETATM 1451  PA  5GP H 187      29.875  44.488  69.823  1.00 19.62
HETATM 1452  O1A 5GP H 187      28.526  44.888  69.143  1.00 19.86
HETATM 1453  O2A 5GP H 187      30.764  44.617  68.702  1.00 23.42
HETATM 1454  O3A 5GP H 187      30.319  45.004  71.073  1.00 20.20
HETATM 1455  O5* 5GP H 187      29.683  43.016  70.027  1.00 20.32
HETATM 1456  C5* 5GP H 187      30.740  42.297  70.837  1.00 21.47
HETATM 1457  C4* 5GP H 187      30.677  40.747  70.770  1.00 21.56
HETATM 1458  O4* 5GP H 187      29.608  40.160  71.599  1.00 20.50
HETATM 1459  C3* 5GP H 187      30.547  40.121  69.352  1.00 20.18
HETATM 1460  O3* 5GP H 187      31.228  38.864  69.416  1.00 23.65
HETATM 1461  C2* 5GP H 187      29.031  39.871  69.248  1.00 18.78
HETATM 1462  O2* 5GP H 187      28.685  38.690  68.496  1.00 20.45
HETATM 1463  C1* 5GP H 187      28.634  39.641  70.688  1.00 17.09
HETATM 1464  N9  5GP H 187      27.238  39.525  71.076  1.00 15.35
HETATM 1465  C8  5GP H 187      26.330  40.535  70.852  1.00 12.57
HETATM 1466  N7' 5GP H 187      25.175  40.314  71.417  1.00 12.88
HETATM 1467  C5  5GP H 187      25.278  39.082  72.070  1.00 10.75
HETATM 1468  C6  5GP H 187      24.326  38.354  72.827  1.00  9.77
HETATM 1469  O6  5GP H 187      23.169  38.678  73.029  1.00  8.66
HETATM 1470  N1' 5GP H 187      24.836  37.190  73.270  1.00  9.67
HETATM 1471  C2  5GP H 187      26.075  36.701  73.001  1.00  9.84
HETATM 1472  N2  5GP H 187      26.361  35.490  73.520  1.00  9.77
HETATM 1473  N3  5GP H 187      27.005  37.353  72.310  1.00 10.31
HETATM 1474  C4  5GP H 187      26.583  38.559  71.844  1.00 12.50
ENDMDL
MODEL        3
HETATM 1475  S   SO4 S 188      31.424  42.923  60.396  1.00 55.69           S4+
HETATM 1476  O1  SO4 S 188      31.631  41.513  60.336  1.00 59.84           O1-
HETATM 1477  O2  SO4 S 188      32.533  43.699  59.932  1.00 49.98           O1-
HETATM 1478  O3  SO4 S 188      31.128  43.217  61.738  1.00 59.44           O1-
HETATM 1479  O4  SO4 S 188      30.353  43.201  59.539  1.00 60.54           O1-
HETATM 1480  O   HOH W 200      29.478  23.354  61.364  1.00  8.67      WATE
ENDMDL
END
""".splitlines()
  stage_1 = pdb.interpretation.stage_1(raw_records=pdb_file)
  sel_cache = stage_1.selection_cache()
  isel = sel_cache.iselection
  assert isel("all").size() == sel_cache.n_seq
  assert isel("none").size() == 0
  assert isel("not all").size() == 0
  assert isel("not none").size() == sel_cache.n_seq
  assert list(isel(r"name c?\*")) == [45,46,48,50,52]
  assert list(isel(r"name 'C?\*'")) == []
  assert list(isel(r"name ' C?\*'")) == [45,46,48,50,52]
  assert list(isel(r"name ' c?\*'")) == []
  assert list(isel(r"name n?'")) == [55, 59]
  assert list(isel(r"altloc a and name n")) == [24]
  assert list(isel(r"altloc b and name n")) == [32]
  assert list(isel(r"altloc ' ' and name n")) == [0,6,17]
  assert list(isel(r"altid ' ' and name n")) == [0,6,17]
  assert list(isel(r"resname hoh")) == [69]
  assert list(isel(r"resname SO4")) == [64,65,66,67,68]
  assert list(isel(r"resname so4")) == [64,65,66,67,68]
  assert list(isel(r"resname So4")) == []
  assert list(isel(r"resname S?4")) == [64,65,66,67,68]
  assert list(isel(r"resname pro and name cg")) == [22]
  assert list(isel(r"resname pro and (name cg or name ca)")) == [18,22]
  assert list(isel(r"not resname pro and (name cg or name ca)")
              ) == [1,7,11,25,33]
  assert list(isel(r"chain h and name o*")) == [41,42,43,44,47,49,51,58]
  assert list(isel(r"(chain h or chain s) and name o[2-46]")) == [58,66,67,68]
  assert list(isel(r"resseq 188")) == [64,65,66,67,68]
  assert list(isel(r"resid 188")) == [64,65,66,67,68]
  assert list(isel(r"resseq 188:")) == [64,65,66,67,68,69]
  assert list(isel(r"resseq 188-")) == [64,65,66,67,68,69]
  assert list(isel(r"resseq :2 and name n*")) == [0,6,13,15,16]
  assert list(isel(r"resseq -2 and name n*")) == [0,6,13,15,16]
  assert list(isel(r"resseq 2:4 and name cb")) == [10,21,28,36]
  assert list(isel(r"resseq 2-4 and name cb")) == [10,21,28,36]
  assert list(isel(r"resid 2-4 and name cb")) == [10,21,28,36]
  assert list(isel(r"model 1 and name cb")) == [4,10,21,28,36]
  assert list(isel(r"model 2-3 and name o1*")) == [41,65]
  assert list(isel(r"icode j and name c?")) == [18,21,22,23]
  assert list(isel(r"segid wate")) == [69]
  assert list(isel(r"element o")) == [65,66,67,68]
  assert list(isel(r"charge 4+")) == [64]

def exercise_xray_structure(anisotropic_flag, verbose=0):
  structure = random_structure.xray_structure(
    space_group_info=sgtbx.space_group_info("P 31"),
    elements=["N","C","C","O","Si"]*2,
    volume_per_atom=500,
    min_distance=2.,
    general_positions_only=False,
    random_u_iso=True,
    anisotropic_flag=anisotropic_flag)
  f_abs = abs(structure.structure_factors(
    anomalous_flag=False, d_min=2, algorithm="direct").f_calc())
  for res_name in (None, "res"):
    for fractional_coordinates in (False, True):
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
  exercise_scalen()
  exercise_selection()
  for anisotropic_flag in (False, True):
    exercise_xray_structure(anisotropic_flag, verbose=verbose)
  print "OK"

if (__name__ == "__main__"):
  run()
