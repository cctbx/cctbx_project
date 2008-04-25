from iotbx import pdb
import iotbx.pdb.remark_290_interpretation
import iotbx.pdb.atom
import iotbx.pdb.parser
from cctbx import crystal
from cctbx import sgtbx
from cctbx.development import random_structure
from cctbx.array_family import flex
import scitbx.math
from libtbx.test_utils import approx_equal, show_diff
import libtbx.load_env
from cStringIO import StringIO
import sys, os

def exercise_combine_unique_pdb_files():
  for file_name,s in [("tmp1", "1"),
                      ("tmp2", "        2"),
                      ("tmp3", "1\t"),
                      ("tmp4", " \t2"),
                      ("tmp5", "1  ")]:
    open(file_name, "w").write(s)
  for file_names in [[], ["tmp1"], ["tmp1", "tmp2"]]:
    c = pdb.combine_unique_pdb_files(file_names=file_names)
    assert len(c.file_name_registry) == len(file_names)
    assert len(c.md5_registry) == len(file_names)
    assert len(c.unique_file_names) == len(file_names)
    assert len(c.raw_records) == len(file_names)
    s = StringIO()
    c.report_non_unique(out=s)
    assert len(s.getvalue()) == 0
  c = pdb.combine_unique_pdb_files(file_names=["tmp1", "tmp1"])
  assert len(c.file_name_registry) == 1
  assert len(c.md5_registry) == 1
  assert len(c.unique_file_names) == 1
  assert len(c.raw_records) == 1
  s = StringIO()
  c.report_non_unique(out=s)
  assert not show_diff(s.getvalue(), """\
INFO: PDB file name appears 2 times: "tmp1"
  1 repeated file name ignored.

""")
  c = pdb.combine_unique_pdb_files(file_names=["tmp1", "tmp1", "tmp2", "tmp1"])
  assert len(c.file_name_registry) == 2
  assert len(c.md5_registry) == 2
  assert len(c.unique_file_names) == 2
  assert len(c.raw_records) == 2
  s = StringIO()
  c.report_non_unique(out=s, prefix="^")
  assert not show_diff(s.getvalue(), """\
^INFO: PDB file name appears 3 times: "tmp1"
^  2 repeated file names ignored.
^
""")
  c = pdb.combine_unique_pdb_files(file_names=["tmp1", "tmp2", "tmp3"])
  assert len(c.file_name_registry) == 3
  assert len(c.md5_registry) == 2
  assert len(c.unique_file_names) == 2
  assert len(c.raw_records) == 2
  s = StringIO()
  c.report_non_unique(out=s)
  assert not show_diff(s.getvalue(), """\
INFO: PDB files with identical content:
  "tmp1"
  "tmp3"
1 file with repeated content ignored.

""")
  c = pdb.combine_unique_pdb_files(file_names=["tmp1", "tmp2", "tmp3", "tmp5"])
  assert len(c.file_name_registry) == 4
  assert len(c.md5_registry) == 2
  assert len(c.unique_file_names) == 2
  assert len(c.raw_records) == 2
  s = StringIO()
  c.report_non_unique(out=s, prefix=": ")
  assert not show_diff(s.getvalue(), """\
: INFO: PDB files with identical content:
:   "tmp1"
:   "tmp3"
:   "tmp5"
: 2 files with repeated content ignored.
:
""")
  c = pdb.combine_unique_pdb_files(file_names=[
    "tmp1", "tmp2", "tmp3", "tmp4", "tmp5", "tmp4", "tmp5"])
  assert len(c.file_name_registry) == 5
  assert len(c.md5_registry) == 2
  assert len(c.unique_file_names) == 2
  assert len(c.raw_records) == 2
  s = StringIO()
  c.report_non_unique(out=s)
  assert not show_diff(s.getvalue(), """\
INFO: PDB file name appears 2 times: "tmp4"
INFO: PDB file name appears 2 times: "tmp5"
  2 repeated file names ignored.
INFO: PDB files with identical content:
  "tmp2"
  "tmp4"
INFO: PDB files with identical content:
  "tmp1"
  "tmp3"
  "tmp5"
3 files with repeated content ignored.

""")

def exercise_pdb_codes_fragment_files():
  all_codes = {} # FUTURE: set
  first_codes = []
  lines = pdb.pdb_codes_fragment_files.splitlines()
  for line in lines:
    codes = line.split()
    for code in codes:
      assert len(code) == 4
      assert not code in all_codes
      all_codes[code] = None
    assert len(codes) >= 2
    assert sorted(codes) == codes
    first_codes.append(codes[0])
  assert sorted(first_codes) == first_codes

def exercise_format_records():
  crystal_symmetry = crystal.symmetry(
    unit_cell=(10,10,13,90,90,120),
    space_group_symbol="R 3").primitive_setting()
  assert iotbx.pdb.format_cryst1_record(crystal_symmetry=crystal_symmetry) \
    == "CRYST1    7.219    7.219    7.219  87.68  87.68  87.68 R 3 :R"
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
  #
  f = iotbx.pdb.format_cryst1_and_scale_records
  assert not show_diff(f(), """\
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1
SCALE1      1.000000  0.000000  0.000000        0.00000
SCALE2      0.000000  1.000000  0.000000        0.00000
SCALE3      0.000000  0.000000  1.000000        0.00000""")
  assert not show_diff(f(crystal_symmetry=crystal_symmetry), """\
CRYST1    7.219    7.219    7.219  87.68  87.68  87.68 R 3 :R
SCALE1      0.138527 -0.005617 -0.005402        0.00000
SCALE2      0.000000  0.138641 -0.005402        0.00000
SCALE3      0.000000  0.000000  0.138746        0.00000""")
  for s in [f(crystal_symmetry=crystal_symmetry.unit_cell()),
            f(crystal_symmetry=crystal_symmetry.unit_cell().parameters()),
            f(scale_fractionalization_matrix=crystal_symmetry.unit_cell()
              .fractionalization_matrix())]:
    assert not show_diff(s, """\
CRYST1    7.219    7.219    7.219  87.68  87.68  87.68 P 1
SCALE1      0.138527 -0.005617 -0.005402        0.00000
SCALE2      0.000000  0.138641 -0.005402        0.00000
SCALE3      0.000000  0.000000  0.138746        0.00000""")
  assert not show_diff(f(cryst1_z=3), """\
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           3
SCALE1      1.000000  0.000000  0.000000        0.00000
SCALE2      0.000000  1.000000  0.000000        0.00000
SCALE3      0.000000  0.000000  1.000000        0.00000""")
  assert not show_diff(f(write_scale_records=False), """\
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1""")
  assert not show_diff(f(scale_u=(1,2,3)), """\
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1
SCALE1      1.000000  0.000000  0.000000        1.00000
SCALE2      0.000000  1.000000  0.000000        2.00000
SCALE3      0.000000  0.000000  1.000000        3.00000""")
  #
  assert iotbx.pdb.format_atom_record() \
    == "ATOM      0  C   DUM     1       0.000   0.000   0.000  1.00  0.00"
  assert iotbx.pdb.format_anisou_record() \
    == "ANISOU    0  C   DUM     1        0      0      0      0      0      0"
  assert iotbx.pdb.format_ter_record() \
    == "TER       0      DUM     1"
  assert iotbx.pdb.format_atom_record(serial=-1, resSeq=-999) \
    == "ATOM     -1  C   DUM  -999       0.000   0.000   0.000  1.00  0.00"
  assert iotbx.pdb.format_anisou_record(serial=100000, resSeq=10000) \
    == "ANISOUA0000  C   DUM  A000        0      0      0      0      0      0"
  assert iotbx.pdb.format_ter_record(serial=100001, resSeq=10001) \
    == "TER   A0001      DUM  A001"
  assert iotbx.pdb.format_atom_record(serial="xyzab", resSeq="TUVW") \
    == "ATOM  xyzab  C   DUM  TUVW       0.000   0.000   0.000  1.00  0.00"

def exercise_parser():
  for i,raw_record in enumerate("""\
LINK         O1  DDA     1                 C3  DDL     2
LINK        MN    MN   391                 OE2 GLU   217            2565
LINK         NZ  LYS A 680        1.260    C4A PLP D   1                LYS-PLP
""".splitlines()):
    r = iotbx.pdb.parser.pdb_record(raw_record=raw_record)
    _1 = [r.name1,r.altLoc1,r.resName1,r.chainID1,r.resSeq1,r.iCode1,r.sym1]
    _2 = [r.name2,r.altLoc2,r.resName2,r.chainID2,r.resSeq2,r.iCode2,r.sym2]
    if (i == 0):
      assert _1 == [' O1 ', ' ', 'DDA', ' ', '   1', ' ', '      ']
      assert _2 == [' C3 ', ' ', 'DDL', ' ', '   2', ' ', '      ']
      assert r.distance is None
      assert r.margin == '        '
    elif (i == 1):
      assert _1 == ['MN  ', ' ', ' MN', ' ', ' 391', ' ', '      ']
      assert _2 == [' OE2', ' ', 'GLU', ' ', ' 217', ' ', '  2565']
      assert r.distance is None
      assert r.margin == '        '
    else:
      assert _1 == [' NZ ', ' ', 'LYS', 'A', ' 680', ' ', '      ']
      assert _2 == [' C4A', ' ', 'PLP', 'D', '   1', ' ', '      ']
      assert r.distance == float(1.260)
      assert r.margin == 'LYS-PLP '

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
  assert str(attr) == " N  ,A,CYS,\",   6, ,1ETN,None"
  s = StringIO()
  attr.show(f=s, prefix="  ")
  assert not show_diff(s.getvalue(), """\
  record name: ATOM
  name:        " N  "
  altLoc:      "A"
  resName:     "CYS"
  chainID:     "\\""
  resSeq:      "   6"
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
""")
  #
  atom_records = """\
ATOM      1  CA  CYS A   6       0.000   0.000   0.000  1.00  0.00
ATOM      2  CA  CYSB    6       0.000   0.000   0.000  1.00  0.00
ATOM      3  CA  CYSAB   7       0.000   0.000   0.000  1.00  0.00
""".splitlines()
  expected_pdb_format = iter("""\
" CA  CYS A   6 "
" CA  CYSB    6 "
" CA  CYSAB   7 "
""".splitlines())
  for record,chainID in zip(atom_records, ["A", "B ", "AB"]):
    attr = pdb.atom.attributes()
    assert attr.pdb_format() == '"               "'
    attr.set_from_ATOM_record(pdb.parser.pdb_record(raw_record=record))
    assert attr.chainID == chainID
    assert attr.pdb_format() == expected_pdb_format.next()

def exercise_residue_name_plus_atom_names_interpreter():
  rnpani = iotbx.pdb.residue_name_plus_atom_names_interpreter
  i = rnpani(residue_name="", atom_names=[])
  assert i.work_residue_name is None
  assert i.atom_name_interpretation is None
  i = rnpani(residue_name="TD", atom_names=["X"])
  assert i.work_residue_name is None
  assert i.atom_name_interpretation is None
  i = rnpani(
    residue_name="thr",
    atom_names=[
      "N", "CA", "C", "O", "CB", "OG1", "CG2",
      "H", "HA", "HB", "HE", "HG1", "1HG2", "2HG2", "3HG2"])
  assert i.work_residue_name == "THR"
  assert i.atom_name_interpretation.unexpected == ["HE"]
  for residue_name in ["c", "dc", "cd", "cyt"]:
    i = rnpani(
      residue_name=residue_name,
      atom_names=[
        "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'",
        "C2'", "C1'", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6",
        "H5'", "H5''", "H4'", "H3'", "H2'", "H2''", "H1'", "H41", "H42",
        "H5", "H6"])
    assert i.work_residue_name == "DC"
    assert i.atom_name_interpretation.unexpected_atom_names() == []

def exercise_format_fasta():
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb",
    test=os.path.isdir)
  if (regression_pdb is None):
    print "Skipping exercise_format_fasta(): input files not available"
    return
  looking_for = ["1ee3_stripped.pdb", "jcm.pdb", "pdb1zff.ent"]
  for node in os.listdir(regression_pdb):
    if (not (node.endswith(".pdb") or node.endswith(".ent"))): continue
    pdb_inp = pdb.input(file_name=os.path.join(regression_pdb, node))
    hierarchy = pdb_inp.construct_hierarchy()
    fasta = []
    for model in hierarchy.models():
      for chain in model.chains():
        for conformer in chain.conformers():
          fasta.append(conformer.format_fasta())
    if (node == "pdb1zff.ent"):
      assert fasta == [
        ['> chain " A" conformer ""', "CCGAATTCGG"]]
      looking_for.remove(node)
    elif (node == "1ee3_stripped.pdb"):
      assert fasta == [
        ['> chain " P" conformer "A"', 'IMEHTV'],
        ['> chain " P" conformer "B"', 'IMEHTV'],
        None,
        None]
      looking_for.remove(node)
    elif (node == "jcm.pdb"):
      assert fasta == [
        ['> chain " A" conformer ""',
         'MSSIFINREYLLPDYIPDELPHREDQIRKIASILAPLYREEKPNNIFIY'
         'GLTGTGKTAVVKFVLSKLHKKFLGKFKHVY',
         'INTRQIDTPYRVLADLLESLDVKVPFTGLSIAELYRRLVKAVRDYGSQV'
         'VIVLDEIDAFVKKYNDDILYKLSRINSISF',
         'IGITNDVKFVDLLDPRVKSSLSEEEIIFPPYNAEELEDILTKRAQMAFK'
         'PGVLPDNVIKLCAALAAREHGDARRALDLL',
         'RVSGEIAERMKDTKVKEEYVYMAKEEIERDRVRDIILTLPFHSKLVLMA'
         'VVSISSEENVVSTTGAVYETYLNICKKLGV',
         'EAVTQRRVSDIINELDMVGILTVVNRGRYGKTKEIGLAVDKNIIVRSLIESDS'],
        ['> chain " B" conformer ""',
         'KNPKVFIDPLSVFKEIPFREDILRDAAIAIRYFVKNEVKFSNLFLGLTG'
         'TGKTFVSKYIFNEIEEVKKEDEEYKDVKQA',
         'YVNCREVGGTPQAVLSSLAGKLAGFSVPKHGINLGEYIDKIKNGTRNIR'
         'AIIYLDEVDTLVKRRGGDIVLYQLLRSDAN',
         'ISVIMISNDINVRDYMEPRVLSSLGPSVIFKPYDAEQLKFILSKYAEYG'
         'LIKGTYDDEILSYIAAISAKEHGDARKAVN',
         'LLFRAAQLASGGGIIRKEHVDKAIVDYEQERLIEAVKALPFHYKLALRS'
         'LIESEDVMSAHKMYTDLCNKFKQKPLSYRR',
         'FSDIISELDMFGIVKIRIINRGRAGGVKKYALVEDKEKVLRALNET'],
        ['> chain " C" conformer ""',
         'TGTAAATTTCCTACGTTTCATCTGAAAATCTAGAGATTTTCAGATGAAACGTAGGAAATTTACATC'],
         None]
      looking_for.remove(node)
  if (len(looking_for) != 0):
    print "WARNING: exercise_format_fasta(): some input files missing:", \
      looking_for

def exercise_xray_structure(use_u_aniso, verbose=0):
  structure = random_structure.xray_structure(
    space_group_info=sgtbx.space_group_info("P 31"),
    elements=["N","C","C","O","Si"]*2,
    volume_per_atom=500,
    min_distance=2.,
    general_positions_only=False,
    random_u_iso=True,
    use_u_aniso=use_u_aniso)
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
      structure_read = iotbx.pdb.input(
        source_info=None,
        lines=flex.std_string(pdb_file.splitlines())).xray_structure_simple(
          fractional_coordinates=fractional_coordinates,
          use_scale_matrix_if_available=False)
      f_read = abs(f_abs.structure_factors_from_scatterers(
        xray_structure=structure_read, algorithm="direct").f_calc())
      regression = flex.linear_regression(f_abs.data(), f_read.data())
      assert regression.is_well_defined()
      if (0 or verbose):
        regression.show_summary()
      assert approx_equal(regression.slope(), 1, eps=1.e-2)
      assert approx_equal(
        regression.y_intercept(), 0, eps=flex.max(f_abs.data())*0.01)

def write_icosahedron():
  for level in xrange(3):
    icosahedron = scitbx.math.icosahedron(level=level)
    scale = 1.5/icosahedron.next_neighbors_distance()
    f = open("icosahedron_%d.pdb"%level, "w")
    for i,site in enumerate(icosahedron.sites*scale):
      print >> f, iotbx.pdb.format_atom_record(serial=i+1, site=site)
    print >> f, "END"
    f.close()

def run():
  verbose = "--verbose" in sys.argv[1:]
  exercise_combine_unique_pdb_files()
  exercise_pdb_codes_fragment_files()
  exercise_format_records()
  exercise_remark_290_interpretation()
  exercise_parser()
  exercise_atom()
  exercise_residue_name_plus_atom_names_interpreter()
  exercise_format_fasta()
  for use_u_aniso in (False, True):
    exercise_xray_structure(use_u_aniso, verbose=verbose)
  write_icosahedron()
  print "OK"

if (__name__ == "__main__"):
  run()
