from iotbx import pdb
import iotbx.pdb.parser
from cctbx.array_family import flex
from libtbx.utils import format_cpu_times
from libtbx.test_utils import show_diff
import libtbx.load_env
import os

def exercise_columns_73_76_evaluator():
  pdb_dir = libtbx.env.find_in_repositories("regression/pdb")
  if (pdb_dir is None):
    print "Skipping exercise_columns_73_76_evaluator():" \
          " input files not available"
    return
  for node in os.listdir(pdb_dir):
    if (not (node.endswith(".pdb") or node.endswith(".ent"))): continue
    file_name = os.path.join(pdb_dir, node)
    lines = flex.split_lines(open(file_name).read())
    py_eval = pdb.parser.columns_73_76_evaluator(raw_records=lines)
    cpp_eval = pdb.columns_73_76_evaluator(lines=lines)
    assert cpp_eval.finding == py_eval.finding
    assert cpp_eval.is_old_style == py_eval.is_old_style

def exercise_line_info_exceptions():
  pdb.input(source_info=None, lines=flex.std_string(["ATOM"]))
  #
  try:
    pdb.input(
      source_info="some.pdb",
      lines=flex.std_string([
        "ATOM   1045  O   HOH  0+30       0.530  42.610  45.267  1.00 33.84"]))
  except RuntimeError, e:
    assert not show_diff(str(e), """\
some.pdb, line 1:
  ATOM   1045  O   HOH  0+30       0.530  42.610  45.267  1.00 33.84
  -----------------------^
  unexpected plus sign.""")
  else: raise RuntimeError("Exception expected.")
  try:
    pdb.input(
      source_info=None,
      lines=flex.std_string([
        "ATOM   1045  O   HOH    30       0.530  42.610  45.267  1.00 33.84",
        "ATOM   1045  O   HOH    3-       0.530  42.610  45.267  1.00 33.84"]))
  except RuntimeError, e:
    assert not show_diff(str(e), """\
input line 2:
  ATOM   1045  O   HOH    3-       0.530  42.610  45.267  1.00 33.84
  -------------------------^
  unexpected minus sign.""")
  else: raise RuntimeError("Exception expected.")
  try:
    pdb.input(
      source_info=None,
      lines=flex.std_string([
        "ATOM   1045  O   HOH    30       0.530  42.610  45.267  1.00 33.84",
        "ATOM   1045  O   HOH    30       0.530  42.610  45.267  1.00 33.84",
        "ATOM   1045  O   HOH  c          0.530  42.610  45.267  1.00 33.84"]))
  except RuntimeError, e:
    assert not show_diff(str(e), """\
input line 3:
  ATOM   1045  O   HOH  c          0.530  42.610  45.267  1.00 33.84
  ----------------------^
  unexpected character.""")
  else: raise RuntimeError("Exception expected.")
  #
  try:
    pdb.input(
      source_info="some.pdb",
      lines=flex.std_string([
        "ATOM   1045  O   HOH    30    x  0.530  42.610  45.267  1.00 33.84"]))
  except RuntimeError, e:
    assert not show_diff(str(e), """\
some.pdb, line 1:
  ATOM   1045  O   HOH    30    x  0.530  42.610  45.267  1.00 33.84
  ------------------------------^
  not a floating-point number.""")
  else: raise RuntimeError("Exception expected.")
  try:
    pdb.input(
      source_info="some.pdb",
      lines=flex.std_string([
        "ATOM   1045  O   HOH    30     x 0.530  42.610  45.267  1.00 33.84"]))
  except RuntimeError, e:
    assert not show_diff(str(e), """\
some.pdb, line 1:
  ATOM   1045  O   HOH    30     x 0.530  42.610  45.267  1.00 33.84
  -------------------------------^
  not a floating-point number.""")
  else: raise RuntimeError("Exception expected.")
  try:
    pdb.input(
      source_info="some.pdb",
      lines=flex.std_string([
        "ATOM   1045  O   HOH    30       0x530  42.610  45.267  1.00 33.84"]))
  except RuntimeError, e:
    assert not show_diff(str(e), """\
some.pdb, line 1:
  ATOM   1045  O   HOH    30       0x530  42.610  45.267  1.00 33.84
  ----------------------------------^
  unexpected character.""")
  else: raise RuntimeError("Exception expected.")

def exercise_pdb_input_process():
  for i_trial in xrange(3):
    pdb_inp = pdb.input(
      source_info=None,
      lines=flex.split_lines(""))
    assert pdb_inp.source_info() == ""
    assert len(pdb_inp.record_type_counts()) == 0
    assert pdb_inp.unknown_section().size() == 0
    assert pdb_inp.title_section().size() == 0
    assert pdb_inp.remark_section().size() == 0
    assert pdb_inp.primary_structure_section().size() == 0
    assert pdb_inp.heterogen_section().size() == 0
    assert pdb_inp.secondary_structure_section().size() == 0
    assert pdb_inp.connectivity_annotation_section().size() == 0
    assert pdb_inp.miscellaneous_features_section().size() == 0
    assert pdb_inp.crystallographic_section().size() == 0
    assert pdb_inp.input_atom_labels_list().size() == 0
    assert pdb_inp.model_numbers().size() == 0
    assert pdb_inp.model_indices().size() == 0
    assert pdb_inp.ter_indices().size() == 0
    assert pdb_inp.chain_indices().size() == 0
    assert pdb_inp.break_indices().size() == 0
    assert pdb_inp.connectivity_section().size() == 0
    assert pdb_inp.bookkeeping_section().size() == 0
    assert pdb_inp.model_numbers_are_unique()
    assert pdb_inp.model_atom_counts().size() == 0
    assert len(pdb_inp.find_duplicate_atom_labels()) == 0
    assert len(pdb_inp.find_false_blank_altloc()) == 0
    pdb_inp = pdb.input(
      source_info="file/name",
      lines=flex.split_lines("""\
HEADER    ISOMERASE                               02-JUL-92   1FKB
ONHOLD    26-JUN-99
OBSLTE     07-DEC-04 1A0Y      1Y4P
TITLE     ATOMIC STRUCTURE OF THE RAPAMYCIN HUMAN IMMUNOPHILIN FKBP-
COMPND    FK506 BINDING PROTEIN (FKBP) COMPLEX WITH IMMUNOSUPPRESSANT
SOURCE    HUMAN (HOMO SAPIENS) RECOMBINANT FORM EXPRESSED IN
KEYWDS    ISOMERASE
EXPDTA    X-RAY DIFFRACTION
AUTHOR    G.D.VAN DUYNE,R.F.STANDAERT,S.L.SCHREIBER,J.C.CLARDY
REVDAT   1   31-OCT-93 1FKB    0
JRNL        AUTH   G.D.VAN DUYNE,R.F.STANDAERT,S.L.SCHREIBER,J.CLARDY
SPRSDE     02-SEP-03 1O58      1J6N
CAVEAT     1B7F    INCORRECT CHIRALITY AT C1* OF U2, CHAIN Q

REMARK   2 RESOLUTION. 1.7  ANGSTROMS.
FTNOTE   1 CIS PEPTIDE: GLY     190  - PHE     191

DBREF  1HTQ A  601   468  SWS    Q10377   GLN1_MYCTU       2    478
SEQRES   1 A  477  THR GLU LYS THR PRO ASP ASP VAL PHE LYS LEU ALA LYS
SEQADV 1KEH ALA A  170  SWS  Q9L5D6    SER   199 ENGINEERED
MODRES 6NSE CYS A  384  CYS  MODIFIED BY CAD

HET    GLC  A 810      12
HETNAM     G6D 6-DEOXY-ALPHA-D-GLUCOSE
HETSYN     G6D QUINOVOSE
FORMUL   2   CA    4(CA1 2+)

HELIX    1   1 GLN A   18  GLY A   34  1                                  17
SHEET    1   A 7 PHE A 257  ALA A 260  0
TURN     1  T1 GLY E   2  THR E   5     BETA, TYPE II

SSBOND  12 CYS B  191    CYS B  220
LINK         N   PRO C  61                 C   GLY A   9            1556
HYDBND       N   GLY A  148                 O   PHE B   41
SLTBRG       N   ILE A  16                 OD2 ASP A 194
CISPEP   1 ALA A  183    PRO A  184          1         0.96

SITE     1 CAB  3 HIS B  57  ASP B 102  SER B 195

CRYST1   45.920   49.790   89.880  90.00  97.34  90.00 P 1 21 1      4
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.021777  0.000000  0.002805        0.00000
SCALE2      0.000000  0.020084  0.000000        0.00000
SCALE3      0.000000  0.000000  0.011218        0.00000
MTRIX1   1  0.739109  0.012922 -0.673462       17.07460    1
MTRIX2   1  0.015672 -0.999875 -0.001986       21.64730    1
MTRIX3   1 -0.673404 -0.009087 -0.739219       44.75290    1
TVECT    1   0.00000   0.00000  20.42000

FOOBAR BAR FOO

MODEL        1
ATOM      1  N   MET A   1       6.215  22.789  24.067  1.00  0.00           N
ATOM      2  CA  MET A   1       6.963  22.789  22.822  1.00  0.00           C
BREAK
HETATM    3  C   MET A   1       7.478  21.387  22.491  1.00  0.00           C
ATOM      4  O   MET A   1       8.406  20.895  23.132  1.00  0.00           O
ENDMDL
MODEL 3
HETATM    9 2H3  MPR     5      16.388   0.289   6.613  1.00  0.08
SIGATM    9 2H3  MPR     5       0.155   0.175   0.155  0.00  0.05
ANISOU    9 2H3  MPR     5      848    848    848      0      0      0
SIGUIJ    9 2H3  MPR     5      510    510    510      0      0      0
TER
ATOM     10  N   CYS     6      14.270   2.464   3.364  1.00  0.07
SIGATM   10  N   CYS     6       0.012   0.012   0.011  0.00  0.00
ANISOU   10  N   CYS     6      788    626    677   -344    621   -232
SIGUIJ   10  N   CYS     6        3     13      4     11      6     13
TER
ENDMDL

CONECT 5332 5333 5334 5335 5336

MASTER       81    0    0    7    3    0    0    645800   20    0   12
END
"""))
    assert pdb_inp.source_info() == "file/name"
    assert pdb_inp.record_type_counts() == {
      "KEYWDS": 1, "SEQRES": 1, "LINK  ": 1, "ORIGX1": 1, "SITE  ": 1,
      "FTNOTE": 1, "HETSYN": 1, "SIGATM": 2, "MTRIX2": 1, "MTRIX3": 1,
      "HELIX ": 1, "MTRIX1": 1, "END   ": 1, "ANISOU": 2, "TITLE ": 1,
      "SLTBRG": 1, "REMARK": 1, "TURN  ": 1, "SCALE1": 1, "SCALE2": 1,
      "AUTHOR": 1, "CRYST1": 1, "SIGUIJ": 2, "CISPEP": 1, "ATOM  ": 4,
      "ENDMDL": 2, "ORIGX2": 1, "MODRES": 1, "SOURCE": 1, "FORMUL": 1,
      "MASTER": 1, "CAVEAT": 1, "HET   ": 1, "COMPND": 1, "MODEL ": 2,
      "REVDAT": 1, "SSBOND": 1, "OBSLTE": 1, "CONECT": 1, "JRNL  ": 1,
      "SPRSDE": 1, "      ":11, "FOOBAR": 1, "HETNAM": 1, "HEADER": 1,
      "ORIGX3": 1, "BREAK ": 1, "ONHOLD": 1, "SHEET ": 1, "TVECT ": 1,
      "HYDBND": 1, "TER   ": 2, "DBREF ": 1, "EXPDTA": 1, "SCALE3": 1,
      "HETATM": 2, "SEQADV": 1}
    assert list(pdb_inp.unknown_section()) == ["FOOBAR BAR FOO"]
    assert not show_diff("\n".join(pdb_inp.title_section()), """\
HEADER    ISOMERASE                               02-JUL-92   1FKB
ONHOLD    26-JUN-99
OBSLTE     07-DEC-04 1A0Y      1Y4P
TITLE     ATOMIC STRUCTURE OF THE RAPAMYCIN HUMAN IMMUNOPHILIN FKBP-
COMPND    FK506 BINDING PROTEIN (FKBP) COMPLEX WITH IMMUNOSUPPRESSANT
SOURCE    HUMAN (HOMO SAPIENS) RECOMBINANT FORM EXPRESSED IN
KEYWDS    ISOMERASE
EXPDTA    X-RAY DIFFRACTION
AUTHOR    G.D.VAN DUYNE,R.F.STANDAERT,S.L.SCHREIBER,J.C.CLARDY
REVDAT   1   31-OCT-93 1FKB    0
JRNL        AUTH   G.D.VAN DUYNE,R.F.STANDAERT,S.L.SCHREIBER,J.CLARDY
SPRSDE     02-SEP-03 1O58      1J6N
CAVEAT     1B7F    INCORRECT CHIRALITY AT C1* OF U2, CHAIN Q""")
    assert not show_diff("\n".join(pdb_inp.remark_section()), """\
REMARK   2 RESOLUTION. 1.7  ANGSTROMS.
FTNOTE   1 CIS PEPTIDE: GLY     190  - PHE     191""")
    assert not show_diff("\n".join(pdb_inp.primary_structure_section()), """\
DBREF  1HTQ A  601   468  SWS    Q10377   GLN1_MYCTU       2    478
SEQRES   1 A  477  THR GLU LYS THR PRO ASP ASP VAL PHE LYS LEU ALA LYS
SEQADV 1KEH ALA A  170  SWS  Q9L5D6    SER   199 ENGINEERED
MODRES 6NSE CYS A  384  CYS  MODIFIED BY CAD""")
    assert not show_diff("\n".join(pdb_inp.heterogen_section()), """\
HET    GLC  A 810      12
HETNAM     G6D 6-DEOXY-ALPHA-D-GLUCOSE
HETSYN     G6D QUINOVOSE
FORMUL   2   CA    4(CA1 2+)""")
    assert not show_diff("\n".join(pdb_inp.secondary_structure_section()), """\
HELIX    1   1 GLN A   18  GLY A   34  1                                  17
SHEET    1   A 7 PHE A 257  ALA A 260  0
TURN     1  T1 GLY E   2  THR E   5     BETA, TYPE II""")
    assert not show_diff(
      "\n".join(pdb_inp.connectivity_annotation_section()), """\
SSBOND  12 CYS B  191    CYS B  220
LINK         N   PRO C  61                 C   GLY A   9            1556
HYDBND       N   GLY A  148                 O   PHE B   41
SLTBRG       N   ILE A  16                 OD2 ASP A 194
CISPEP   1 ALA A  183    PRO A  184          1         0.96""")
    assert not show_diff(
      "\n".join(pdb_inp.miscellaneous_features_section()), """\
SITE     1 CAB  3 HIS B  57  ASP B 102  SER B 195""")
    assert not show_diff("\n".join(pdb_inp.crystallographic_section()), """\
CRYST1   45.920   49.790   89.880  90.00  97.34  90.00 P 1 21 1      4
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.021777  0.000000  0.002805        0.00000
SCALE2      0.000000  0.020084  0.000000        0.00000
SCALE3      0.000000  0.000000  0.011218        0.00000
MTRIX1   1  0.739109  0.012922 -0.673462       17.07460    1
MTRIX2   1  0.015672 -0.999875 -0.001986       21.64730    1
MTRIX3   1 -0.673404 -0.009087 -0.739219       44.75290    1
TVECT    1   0.00000   0.00000  20.42000""")
    assert pdb_inp.input_atom_labels_list().size() == 6
    assert list(pdb_inp.model_numbers()) == [1,3]
    assert list(pdb_inp.model_indices()) == [4,6]
    assert list(pdb_inp.ter_indices()) == [5,6]
    assert [list(v) for v in pdb_inp.chain_indices()] == [[4],[5,6]]
    assert list(pdb_inp.break_indices()) == [2]
    assert not show_diff("\n".join(pdb_inp.connectivity_section()), """\
CONECT 5332 5333 5334 5335 5336""")
    assert not show_diff("\n".join(pdb_inp.bookkeeping_section()), """\
MASTER       81    0    0    7    3    0    0    645800   20    0   12
END""")
    assert pdb_inp.name_selection_cache().keys() \
        == [" C  ", " CA ", " N  ", " O  ", "2H3 "]
    assert [list(v) for v in pdb_inp.name_selection_cache().values()] \
        == [[2], [1], [0,5], [3], [4]]
    assert pdb_inp.altloc_selection_cache().keys() == [" "]
    assert [list(v) for v in pdb_inp.altloc_selection_cache().values()] \
        == [[0,1,2,3,4,5]]
    assert pdb_inp.resname_selection_cache().keys() == ["CYS ", "MET ", "MPR "]
    assert [list(v) for v in pdb_inp.resname_selection_cache().values()] \
        == [[5], [0,1,2,3], [4]]
    assert pdb_inp.chain_selection_cache().keys() == [" ", "A"]
    assert [list(v) for v in pdb_inp.chain_selection_cache().values()] \
        == [[4,5], [0,1,2,3]]
    for resseq,i_seqs in [(1,[0,1,2,3]),(5,[4]),(6,[5])]:
      assert list(pdb_inp.resseq_selection_cache()[resseq+999]) == i_seqs
    for i,i_seqs in enumerate(pdb_inp.resseq_selection_cache()):
      if (i not in [j+999 for j in [1,5,6]]):
        assert i_seqs.size() == 0
    assert pdb_inp.icode_selection_cache().keys() == [" "]
    assert [list(v) for v in pdb_inp.icode_selection_cache().values()] \
        == [[0,1,2,3,4,5]]
    assert pdb_inp.segid_selection_cache().keys() == ["    "]
    assert [list(v) for v in pdb_inp.segid_selection_cache().values()] \
        == [[0,1,2,3,4,5]]
    assert pdb_inp.model_numbers_are_unique()
    assert list(pdb_inp.model_atom_counts()) == [4,2]
    assert len(pdb_inp.find_duplicate_atom_labels()) == 0
    assert len(pdb_inp.find_false_blank_altloc()) == 0
    print pdb_inp.conformer_selections()
  #
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM      1  N   MET A   1       6.215  22.789  24.067  1.00  0.00           N
ATOM      2  N   MET A   1       2.615  27.289  20.467  1.00  0.00           O
"""))
  dup = pdb_inp.find_duplicate_atom_labels()
  assert dup.size() == 1
  assert list(dup[0]) == [0,1]
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
MODEL 1
ATOM      1  N   MET A   1       6.215  22.789  24.067  1.00  0.00           N
ATOM      2  N   MET A   1       2.615  27.289  20.467  1.00  0.00           O
ENDMDL
MODEL 2
ATOM      1  N   MET A   1       6.215  22.789  24.067  1.00  0.00           N
ATOM      2  N   MET A   1       2.615  27.289  20.467  1.00  0.00           O
ENDMDL
"""))
  dup = pdb_inp.find_duplicate_atom_labels()
  assert dup.size() == 2
  assert list(dup[0]) == [0,1]
  assert list(dup[1]) == [2,3]
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM      1  N   MET A   1       6.215  22.789  24.067  1.00  0.00           N
ATOM      2  N   MET A   2       2.615  27.289  20.467  1.00  0.00           O
ATOM      3  N   MET A   1       2.615  27.289  20.467  1.00  0.00           O
ATOM      4  N   MET A   1       2.615  27.289  20.467  1.00  0.00           O
ATOM      5  C   MET A   1       6.215  22.789  24.067  1.00  0.00           N
ATOM      6  C   MET A   2       2.615  27.289  20.467  1.00  0.00           O
ATOM      7  C   MET A   1       2.615  27.289  20.467  1.00  0.00           O
ATOM      8  C   MET A   1       2.615  27.289  20.467  1.00  0.00           O
"""))
  dup = pdb_inp.find_duplicate_atom_labels()
  assert dup.size() == 2
  assert list(dup[0]) == [0,2,3]
  assert list(dup[1]) == [4,6,7]
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM      1  CB  LYS   109      16.113   7.345  47.084  1.00 20.00      A
ATOM      2  CG  LYS   109      17.058   6.315  47.703  1.00 20.00      A
ATOM      3  CB  LYS   109      26.721   1.908  15.275  1.00 20.00      B
ATOM      4  CG  LYS   109      27.664   2.793  16.091  1.00 20.00      B
"""))
  assert pdb_inp.find_duplicate_atom_labels().size() == 0
  expected_pdb_format = """\
" CB  LYS   109 " segid="A   "
" CG  LYS   109 " segid="A   "
" CB  LYS   109 " segid="B   "
" CG  LYS   109 " segid="B   "
""".splitlines()
  for il,ef in zip(pdb_inp.input_atom_labels_list(),expected_pdb_format):
    assert il.pdb_format() == ef
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM  12345qN123AR123C1234Ixyz1234.6781234.6781234.678123.56213.56abcdefS123E1C1
HETATM12345qN123AR123C1234Ixyz1234.6781234.6781234.678123.56213.56abcdefS123E1C1
"""))
  for ial in pdb_inp.input_atom_labels_list():
    assert ial.name() == "N123"
    assert ial.altloc() == "A"
    assert ial.resname() == "R123"
    assert ial.chain() == "C"
    assert ial.resseq == 1234
    assert ial.icode() == "I"
    assert ial.segid() == "S123"
    assert ial.pdb_format() == '"N123AR123C1234I" segid="S123"'
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM  12345qN123AR123C1234Ixyz1234.6781234.6781234.678123.56213.56abcdef    E1C1
HETATM12345qN123AR123C1234Ixyz1234.6781234.6781234.678123.56213.56abcdef    E1C1
"""))
  for ial in pdb_inp.input_atom_labels_list():
    assert ial.name() == "N123"
    assert ial.altloc() == "A"
    assert ial.resname() == "R123"
    assert ial.chain() == "C"
    assert ial.resseq == 1234
    assert ial.icode() == "I"
    assert ial.segid() == "    "
    assert ial.pdb_format() == '"N123AR123C1234I"'
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM
HETATM
"""))
  for ial in pdb_inp.input_atom_labels_list():
    assert ial.name() == "    "
    assert ial.altloc() == " "
    assert ial.resname() == "    "
    assert ial.chain() == " "
    assert ial.resseq == 0
    assert ial.icode() == " "
    assert ial.segid() == "    "
    assert ial.pdb_format() == '"             0 "'
  #
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
"""))
  assert list(pdb_inp.model_indices()) == []
  assert list(pdb_inp.chain_indices()) == []
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM
"""))
  assert list(pdb_inp.model_indices()) == [1]
  assert [list(v) for v in pdb_inp.chain_indices()] == [[1]]
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
MODEL        1
ENDMDL
"""))
  assert list(pdb_inp.model_indices()) == [0]
  assert [list(v) for v in pdb_inp.chain_indices()] == [[]]
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
MODEL        1
ATOM
ENDMDL
"""))
  assert list(pdb_inp.model_indices()) == [1]
  assert [list(v) for v in pdb_inp.chain_indices()] == [[1]]
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
MODEL        1
ENDMDL
MODEL        2
ENDMDL
"""))
  assert list(pdb_inp.model_indices()) == [0,0]
  assert [list(v) for v in pdb_inp.chain_indices()] == [[],[]]
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
MODEL        1
ENDMDL
MODEL        2
ATOM
ENDMDL
"""))
  assert list(pdb_inp.model_indices()) == [0,1]
  assert [list(v) for v in pdb_inp.chain_indices()] == [[],[1]]
  try:
    pdb.input(
      source_info=None,
      lines=flex.split_lines("""\
MODEL        1
ENDMDL
ATOM
"""))
  except RuntimeError, e:
    assert not show_diff(str(e), """\
input line 3:
  ATOM
  ^
  ATOM or HETATM record is outside MODEL/ENDMDL block.""")
  else: raise RuntimeError("Exception expected.")
  try:
    pdb.input(
      source_info=None,
      lines=flex.split_lines("""\
MODEL        1
MODEL        2
"""))
  except RuntimeError, e:
    assert not show_diff(str(e), """\
input line 2:
  MODEL        2
  ^
  Missing ENDMDL for previous MODEL record.""")
  else: raise RuntimeError("Exception expected.")
  try:
    pdb.input(
      source_info=None,
      lines=flex.split_lines("""\
ATOM
MODEL        1
"""))
  except RuntimeError, e:
    assert not show_diff(str(e), """\
input line 2:
  MODEL        1
  ^
  MODEL record must appear before any ATOM or HETATM records.""")
  else: raise RuntimeError("Exception expected.")
  try:
    pdb.input(
      source_info=None,
      lines=flex.split_lines("""\
ATOM
ENDMDL
"""))
  except RuntimeError, e:
    assert not show_diff(str(e), """\
input line 2:
  ENDMDL
  ^
  No matching MODEL record.""")
  else: raise RuntimeError("Exception expected.")
  #
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
MODEL        1
ATOM                 C
ATOM                 C
ATOM                 D
ATOM                 E
ATOM                 E
ENDMDL
MODEL        2
ATOM                 C
ATOM                 C
ATOM                 D
ATOM                 D
ATOM                 E
ENDMDL
MODEL        3
ATOM                                                                    C
ATOM                                                                    D
ATOM                                                                    D
ATOM                 E                                                  X
ATOM                 E
ENDMDL
"""))
  assert list(pdb_inp.model_indices()) == [5,10,15]
  assert [list(v) for v in pdb_inp.chain_indices()] \
      == [[2,3,5],[7,9,10],[11,13,15]]

def exercise():
  exercise_columns_73_76_evaluator()
  exercise_line_info_exceptions()
  exercise_pdb_input_process()
  print format_cpu_times()

if (__name__ == "__main__"):
  exercise()
