from iotbx import pdb
from iotbx.pdb import hybrid_36
from cctbx import crystal
from cctbx.array_family import flex
from libtbx.utils import Sorry, hashlib_md5, \
  user_plus_sys_time, format_cpu_times
from libtbx.test_utils import Exception_expected, approx_equal, show_diff
import libtbx.load_env
from cStringIO import StringIO
import cPickle
import pickle
import sys, os

def exercise_hybrid_36():
  hybrid_36.exercise(hy36enc=pdb.hy36encode, hy36dec=pdb.hy36decode)
  for width,s in [(3,"AAA"), (6,"zzzzzz")]:
    try: pdb.hy36encode(width=width, value=0)
    except RuntimeError, e:
      assert str(e) == "unsupported width."
    else: raise Exception_expected
    try: pdb.hy36decode(width=width, s=s)
    except RuntimeError, e:
      assert str(e) == "unsupported width."
    else: raise Exception_expected
  ups = user_plus_sys_time()
  n_ok = pdb.hy36recode_width_4_all()
  ups = ups.elapsed()
  print "time hy36recode_width_4_all: %.2f s" \
    " (%.3f micro s per encode-decode cycle)" % (ups, 1.e6*ups/max(1,n_ok))
  assert n_ok == 999+10000+2*26*36**3
  #
  assert pdb.resseq_decode(s=1234) == 1234
  assert pdb.resseq_decode(s="A123") == 11371
  assert pdb.resseq_decode(s="1") == 1
  pdb.resseq_encode(value=1) == "   1"
  pdb.resseq_encode(value=11371) == "A123"
  pdb.resseq_encode(value=1234) == "1234"
  #
  try: pdb.resseq_decode(s="18A")
  except ValueError, e:
    assert str(e) == 'invalid residue sequence number: " 18A"'
  else: raise Exception_expected

def exercise_base_256_ordinal():
  o = pdb.utils_base_256_ordinal
  assert o(None) == 48
  assert o("") == 48
  assert o("0") == 48
  assert o(" 0") == 48
  assert o("  0") == 48
  assert o("1") == 49
  assert o("-1") == -49
  assert o("123") == (49*256+50)*256+51
  assert o("-123") == -o("123")
  #
  def po(s):
    result = 0
    s = s.lstrip()
    neg = False
    if (len(s) == 0):
      s = "0"
    elif (s[0] == "-"):
      neg = True
      s = s[1:]
    for c in s:
      result *= 256
      result += ord(c)
    if (neg):
      result *= -1
    return result
  for s in ["780 ", "999 ", "1223 ", "zzzzz"]:
    assert o(s) == po(s)
    assert o("-"+s) == -o(s)
  #
  def o_cmp(a, b): return cmp(o(a), o(b))
  char4s = ["%4s" % i for i in xrange(-999,9999+1)]
  assert sorted(char4s, o_cmp) == char4s
  m = pdb.hy36decode(width=4, s="zzzz")
  e = pdb.hy36encode
  char4s = [e(width=4, value=i) for i in xrange(-999,m+1,51)]
  assert sorted(char4s, o_cmp) == char4s

def exercise_columns_73_76_evaluator(pdb_file_names):
  if (pdb_file_names is None):
    print "Skipping exercise_columns_73_76_evaluator():" \
          " input files not available"
    return
  known_blank = """\
occ_3_bad2.pdb
enk_gm.pdb
t.pdb
phe_a.pdb
f_obs_complex.pdb
phe_h_bad.pdb
one_conf_but_altloc.pdb
""".splitlines()
  known_exactly_one = """\
pdb103l.ent
pdb1etn.ent
pdb118d.ent
pdb161d.ent
pdb139l.ent
pdb1anp.ent
pdb1gky.ent
""".splitlines()
  n_known = [0, 0]
  for file_name in pdb_file_names:
    lines = flex.split_lines(open(file_name).read())
    e = pdb.columns_73_76_evaluator(lines=lines)
    bn = os.path.basename(file_name)
    if (bn in known_blank):
      assert e.finding == "Blank columns 73-76 on ATOM and HETATM records."
      assert not e.is_old_style
      n_known[0] += 1
    elif (bn in known_exactly_one):
      assert e.finding == "Exactly one common label in columns 73-76."
      assert e.is_old_style
      n_known[1] += 1
  assert n_known[0] >= 3
  assert n_known[1] >= 3
  #
  lines = flex.split_lines("""\
HEADER    HYDROLASE(METALLOPROTEINASE)            17-NOV-93   1THL
ATOM      1  N   ILE     1       9.581  51.813  -0.720  1.00 31.90      1THL 158
ATOM      2  CA  ILE     1       8.335  52.235  -0.041  1.00 52.95      1THL 159
ATOM      3  C   ILE     1       7.959  53.741   0.036  1.00 26.88      1THL 160
END
""")
  e = pdb.columns_73_76_evaluator(lines=lines)
  assert e.finding == "Exactly one common label in columns 73-76."
  assert e.is_old_style

def exercise_line_info_exceptions():
  pdb.input(source_info=None, lines=flex.std_string(["ATOM"]))
  #
  try:
    pdb.input(
      source_info="some.pdb",
      lines=flex.split_lines("""\
HETATM    9 2H3  MPR B   5      16.388   0.289   6.613  1.00  0.08
ANISOU    9 2H3  MPR B   5      8+8    848    848      0      0      0
"""))
  except ValueError, e:
    assert not show_diff(str(e), """\
some.pdb, line 2:
  ANISOU    9 2H3  MPR B   5      8+8    848    848      0      0      0
  ---------------------------------^
  unexpected plus sign.""")
  else: raise Exception_expected
  try:
    pdb.input(
      source_info=None,
      lines=flex.split_lines("""\
HETATM    9 2H3  MPR B   5      16.388   0.289   6.613  1.00  0.08
HETATM    9 2H3  MPR B   5      16.388   0.289   6.613  1.00  0.08
ANISOU    9 2H3  MPR B   5      84-    848    848      0      0      0
"""))
  except ValueError, e:
    assert not show_diff(str(e), """\
input line 3:
  ANISOU    9 2H3  MPR B   5      84-    848    848      0      0      0
  ----------------------------------^
  unexpected minus sign.""")
  else: raise Exception_expected
  try:
    pdb.input(
      source_info=None,
      lines=flex.split_lines("""\
HETATM    9 2H3  MPR B   5      16.388   0.289   6.613  1.00  0.08
ANISOU    9 2H3  MPR B   5    c        848    848      0      0      0
"""))
  except ValueError, e:
    assert not show_diff(str(e), """\
input line 2:
  ANISOU    9 2H3  MPR B   5    c        848    848      0      0      0
  ------------------------------^
  unexpected character.""")
  else: raise Exception_expected
  #
  try:
    pdb.input(
      source_info="some.pdb",
      lines=flex.std_string([
        "ATOM   1045  O   HOH    30    x  0.530  42.610  45.267  1.00 33.84"]))
  except ValueError, e:
    assert not show_diff(str(e), """\
some.pdb, line 1:
  ATOM   1045  O   HOH    30    x  0.530  42.610  45.267  1.00 33.84
  ------------------------------^
  not a floating-point number.""")
  else: raise Exception_expected
  try:
    pdb.input(
      source_info="some.pdb",
      lines=flex.std_string([
        "ATOM   1045  O   HOH    30     x 0.530  42.610  45.267  1.00 33.84"]))
  except ValueError, e:
    assert not show_diff(str(e), """\
some.pdb, line 1:
  ATOM   1045  O   HOH    30     x 0.530  42.610  45.267  1.00 33.84
  -------------------------------^
  not a floating-point number.""")
  else: raise Exception_expected
  try:
    pdb.input(
      source_info="some.pdb",
      lines=flex.std_string([
        "HETATM 4160  O   HOH S 272         nan   0.000   0.000  1.00 54.72"]))
  except ValueError, e:
    assert not show_diff(str(e), """\
some.pdb, line 1:
  HETATM 4160  O   HOH S 272         nan   0.000   0.000  1.00 54.72
  -----------------------------------^
  not a floating-point number.""")
  else: raise Exception_expected
  try:
    pdb.input(
      source_info="some.pdb",
      lines=flex.std_string([
        "ATOM   1045  O   HOH    30       0x530  42.610  45.267  1.00 33.84"]))
  except ValueError, e:
    assert not show_diff(str(e), """\
some.pdb, line 1:
  ATOM   1045  O   HOH    30       0x530  42.610  45.267  1.00 33.84
  ----------------------------------^
  unexpected character.""")
  else: raise Exception_expected

pdb_string_all_sections = """\
HEADER    ISOMERASE                               02-JUL-92   1FKB
ONHOLD    26-JUN-99
OBSLTE     07-DEC-04 1A0Y      1Y4P
TITLE     ATOMIC STRUCTURE OF THE RAPAMYCIN HUMAN IMMUNOPHILIN FKBP-
SPLIT      2QNH 1VSP
COMPND    FK506 BINDING PROTEIN (FKBP) COMPLEX WITH IMMUNOSUPPRESSANT
SOURCE    HUMAN (HOMO SAPIENS) RECOMBINANT FORM EXPRESSED IN
KEYWDS    ISOMERASE
EXPDTA    X-RAY DIFFRACTION
NUMMDL    8
MDLTYP    CA ATOMS ONLY, CHAIN B, C, D, E, F, G, H, I, J, K, L, M, N,
MDLTYP   2 O, P, Q, R, S, T, U
AUTHOR    G.D.VAN DUYNE,R.F.STANDAERT,S.L.SCHREIBER,J.C.CLARDY
REVDAT   1   31-OCT-93 1FKB    0
JRNL        AUTH   G.D.VAN DUYNE,R.F.STANDAERT,S.L.SCHREIBER,J.CLARDY
SPRSDE     02-SEP-03 1O58      1J6N
CAVEAT     1B7F    INCORRECT CHIRALITY AT C1* OF U2, CHAIN Q

REMARK   2 RESOLUTION. 1.7  ANGSTROMS.
FTNOTE   1 CIS PEPTIDE: GLY     190  - PHE     191

DBREF  1HTQ A  601   468  SWS    Q10377   GLN1_MYCTU       2    478
DBREF1 1JZX A    1  2880  GB                   15805042
DBREF2 1JZX A     NC_001263                     2587937     2590817
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
HETATM    3  C   MET A   2       7.478  21.387  22.491  1.00  0.00           C
ATOM      4  O   MET A   2       8.406  20.895  23.132  1.00  0.00           O
ENDMDL
MODEL 3
HETATM    9 2H3  MPR B   5      16.388   0.289   6.613  1.00  0.08
SIGATM    9 2H3  MPR B   5       0.155   0.175   0.155  0.00  0.05
ANISOU    9 2H3  MPR B   5      848    848    848      0      0      0
SIGUIJ    9 2H3  MPR B   5      510    510    510      0      0      0
TER
ATOM     10  N   CYSCH   6      14.270   2.464   3.364  1.00  0.07
SIGATM   10  N   CYSCH   6       0.012   0.012   0.011  0.00  0.00
ANISOU   10  N   CYSCH   6      788    626    677   -344    621   -232
SIGUIJ   10  N   CYSCH   6        3     13      4     11      6     13
TER
ENDMDL

CONECT 5332 5333 5334 5335 5336

MASTER       81    0    0    7    3    0    0    645800   20    0   12
END
"""

def exercise_pdb_input():
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
    assert len(pdb_inp.atoms_with_labels()) == 0
    assert pdb_inp.atoms().size() == 0
    assert pdb_inp.model_ids().size() == 0
    assert pdb_inp.model_indices().size() == 0
    assert pdb_inp.ter_indices().size() == 0
    assert pdb_inp.chain_indices().size() == 0
    assert pdb_inp.break_indices().size() == 0
    assert pdb_inp.connectivity_section().size() == 0
    assert pdb_inp.bookkeeping_section().size() == 0
    assert pdb_inp.model_atom_counts().size() == 0
    pdb_inp = pdb.input(
      source_info="file/name",
      lines=pdb_string_all_sections)
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
      "HETATM": 2, "SEQADV": 1, "SPLIT ": 1, "NUMMDL": 1, "MDLTYP": 2,
      "DBREF1": 1, "DBREF2": 1}
    assert list(pdb_inp.unknown_section()) == ["FOOBAR BAR FOO"]
    assert not show_diff("\n".join(pdb_inp.title_section()), """\
HEADER    ISOMERASE                               02-JUL-92   1FKB
ONHOLD    26-JUN-99
OBSLTE     07-DEC-04 1A0Y      1Y4P
TITLE     ATOMIC STRUCTURE OF THE RAPAMYCIN HUMAN IMMUNOPHILIN FKBP-
SPLIT      2QNH 1VSP
COMPND    FK506 BINDING PROTEIN (FKBP) COMPLEX WITH IMMUNOSUPPRESSANT
SOURCE    HUMAN (HOMO SAPIENS) RECOMBINANT FORM EXPRESSED IN
KEYWDS    ISOMERASE
EXPDTA    X-RAY DIFFRACTION
NUMMDL    8
MDLTYP    CA ATOMS ONLY, CHAIN B, C, D, E, F, G, H, I, J, K, L, M, N,
MDLTYP   2 O, P, Q, R, S, T, U
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
DBREF1 1JZX A    1  2880  GB                   15805042
DBREF2 1JZX A     NC_001263                     2587937     2590817
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
    assert len(pdb_inp.atoms_with_labels()) == 6
    assert [atom.serial for atom in pdb_inp.atoms()] \
        == ["    1", "    2", "    3", "    4", "    9", "   10"]
    assert [atom.element for atom in pdb_inp.atoms()] \
        == [" N", " C", " C", " O", "  ", "  "]
    assert list(pdb_inp.model_ids()) == ["   1", "   3"]
    assert list(pdb_inp.model_indices()) == [4,6]
    assert list(pdb_inp.ter_indices()) == [5,6]
    assert [list(v) for v in pdb_inp.chain_indices()] == [[4],[5,6]]
    assert list(pdb_inp.break_indices()) == [2]
    assert not show_diff("\n".join(pdb_inp.connectivity_section()), """\
CONECT 5332 5333 5334 5335 5336""")
    assert not show_diff("\n".join(pdb_inp.bookkeeping_section()), """\
MASTER       81    0    0    7    3    0    0    645800   20    0   12
END""")
    assert list(pdb_inp.model_atom_counts()) == [4,2]
  #
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM      1  CB  LYS   109      16.113   7.345  47.084  1.00 20.00      A
ATOM      2  CG  LYS   109      17.058   6.315  47.703  1.00 20.00      A
ATOM      3  CB  LYS   109      26.721   1.908  15.275  1.00 20.00      B
ATOM      4  CG  LYS   109      27.664   2.793  16.091  1.00 20.00      B
"""))
  expected_id_strs = """\
pdb=" CB  LYS   109 " segid="A   "
pdb=" CG  LYS   109 " segid="A   "
pdb=" CB  LYS   109 " segid="B   "
pdb=" CG  LYS   109 " segid="B   "
""".splitlines()
  for awl,eids in zip(pdb_inp.atoms_with_labels(), expected_id_strs):
    assert not show_diff(awl.id_str(), eids)
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM  12345qN123AR12 C1234Ixyz1234.6781234.6781234.678123.56213.56abcdefS123E1C1
HETATM12345qN123AR12 C1234Ixyz1234.6781234.6781234.678123.56213.56abcdefS123E1C1
"""))
  for awl in pdb_inp.atoms_with_labels():
    assert awl.name == "N123"
    assert awl.altloc == "A"
    assert awl.resname == "R12"
    assert awl.chain_id == "C"
    assert awl.resseq == "1234"
    assert awl.icode == "I"
    assert awl.segid == "S123"
    assert awl.id_str() == 'pdb="N123AR12 C1234I" segid="S123"'
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM  12345qN123AR12 C1234Ixyz1234.6781234.6781234.678123.56213.56abcdef    E1C1
HETATM12345qN123AR12 C1234Ixyz1234.6781234.6781234.678123.56213.56abcdef    E1C1
"""))
  for awl in pdb_inp.atoms_with_labels():
    assert awl.name == "N123"
    assert awl.altloc == "A"
    assert awl.resname == "R12"
    assert awl.chain_id == "C"
    assert awl.resseq == "1234"
    assert awl.icode == "I"
    assert awl.segid == "    "
    assert awl.id_str() == 'pdb="N123AR12 C1234I"'
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM
HETATM
"""))
  for awl in pdb_inp.atoms_with_labels():
    assert awl.name == "    "
    assert awl.altloc == " "
    assert awl.resname == "   "
    assert awl.chain_id == " "
    assert awl.resseq == "    "
    assert awl.icode == " "
    assert awl.segid == "    "
    assert awl.id_str() == 'pdb="               "'
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
  except ValueError, e:
    assert not show_diff(str(e), """\
input line 3:
  ATOM
  ^
  ATOM or HETATM record is outside MODEL/ENDMDL block.""")
  else: raise Exception_expected
  try:
    pdb.input(
      source_info=None,
      lines=flex.split_lines("""\
MODEL        1
MODEL        2
"""))
  except ValueError, e:
    assert not show_diff(str(e), """\
input line 2:
  MODEL        2
  ^
  Missing ENDMDL for previous MODEL record.""")
  else: raise Exception_expected
  try:
    pdb.input(
      source_info=None,
      lines=flex.split_lines("""\
ATOM
MODEL        1
"""))
  except ValueError, e:
    assert not show_diff(str(e), """\
input line 2:
  MODEL        1
  ^
  MODEL record must appear before any ATOM or HETATM records.""")
  else: raise Exception_expected
  try:
    pdb.input(
      source_info=None,
      lines=flex.split_lines("""\
ATOM
ENDMDL
"""))
  except ValueError, e:
    assert not show_diff(str(e), """\
input line 2:
  ENDMDL
  ^
  No matching MODEL record.""")
  else: raise Exception_expected
  #
  for record_name in ["SIGATM", "ANISOU", "SIGUIJ"]:
    try:
      pdb.input(source_info=None, lines=flex.std_string([record_name]))
    except ValueError, e:
      assert not show_diff(str(e), """\
input line 1:
  %s
  ^
  no matching ATOM or HETATM record.""" % record_name)
    else: raise Exception_expected
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
MODEL        4
ATOM                                                                    C
ATOM                                                                    D
ATOM                                                                    C
ATOM                 E                                                  X
ATOM                 E
ENDMDL
"""))
  assert list(pdb_inp.model_indices()) == [5,10,15,20]
  assert [list(v) for v in pdb_inp.chain_indices()] \
      == [[2,3,5],[7,9,10],[11,13,15], [18,20]]
  #
  open("tmp.pdb", "w")
  pdb_inp = pdb.input(file_name="tmp.pdb")
  assert pdb_inp.source_info() == "file tmp.pdb"
  open("tmp.pdb", "w").write("""\
ATOM      1  CA  SER     1       1.212 -12.134   3.757  1.00  0.00
ATOM      2  CA  LEU     2       1.118  -9.777   0.735  1.00  0.00
""")
  pdb_inp = pdb.input(file_name="tmp.pdb")
  try: pdb.input(file_name="")
  except IOError, e:
    assert str(e).startswith('Cannot open file for reading: ""')
  else: raise Exception_expected
  #
  assert "HIS" in pdb.common_residue_names_amino_acid
  assert "GUA" in pdb.common_residue_names_rna_dna
  assert "CD " in pdb.common_residue_names_ccp4_mon_lib_rna_dna
  assert "HOH" in pdb.common_residue_names_water
  assert "SO4" in pdb.common_residue_names_small_molecule
  assert " FE" in pdb.common_residue_names_element
  get_class = pdb.common_residue_names_get_class
  assert get_class(name="ALA") == "common_amino_acid"
  assert get_class(name="  U") == "common_rna_dna"
  assert get_class(name="HOH") == "common_water"
  assert get_class(name="SO4") == "common_small_molecule"
  assert get_class(name="CL ") == "common_element"
  assert get_class(name="ABC") == "other"
  assert get_class(name="CD ") == "common_element"
  assert get_class(name="CD ", consider_ccp4_mon_lib_rna_dna=True) \
    == "ccp4_mon_lib_rna_dna"
  #
  assert pdb.rna_dna_reference_residue_name(common_name="ALA") is None
  for common_names in [pdb.common_residue_names_rna_dna,
                       pdb.common_residue_names_ccp4_mon_lib_rna_dna]:
    for n in common_names:
      r = pdb.rna_dna_reference_residue_name(common_name=n)
      assert r is not None
      assert pdb.rna_dna_reference_residue_name(
        common_name=" "+n.lower()+" ") == r
  #
  for line in """\
CRYST1   61.410   54.829   43.543  90.00  90.00  90.00 P 21 21 21    8
REMARK sg= P2(1)2(1)2(1) a= 61.410 b= 54.829 c= 43.543 alpha= 90 beta= 90 gamma= 90
""".splitlines():
    pdb_inp = pdb.input(source_info=None, lines=flex.std_string([line]))
    cs = pdb_inp.crystal_symmetry()
    assert str(cs.unit_cell()) == "(61.41, 54.829, 43.543, 90, 90, 90)"
    assert str(cs.space_group_info()) == "P 21 21 21"
    sps = pdb_inp.special_position_settings()
    assert sps.is_similar_symmetry(cs)
    assert approx_equal(sps.min_distance_sym_equiv(), 0.5)
    for weak_symmetry in [False, True]:
      cs = pdb_inp.crystal_symmetry(
        crystal_symmetry=crystal.symmetry(
          unit_cell=(10,10,10,90,90,90)),
          weak_symmetry=weak_symmetry)
      if (weak_symmetry):
        assert str(cs.unit_cell()) == "(61.41, 54.829, 43.543, 90, 90, 90)"
      else:
        assert str(cs.unit_cell()) == "(10, 10, 10, 90, 90, 90)"
      assert str(cs.space_group_info()) == "P 21 21 21"
      sps = pdb_inp.special_position_settings(
        special_position_settings=cs.special_position_settings(
          min_distance_sym_equiv=3),
        weak_symmetry=weak_symmetry)
      assert sps.is_similar_symmetry(cs)
      assert approx_equal(sps.min_distance_sym_equiv(), 3)
  #
  assert pdb_inp.extract_header_year() is None
  assert pdb_inp.extract_remark_iii_records(iii=2) == []
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
HEADER                                            02-JUL-92
REMARK   2 RESOLUTION. 1.7  ANGSTROMS.
"""))
  assert pdb_inp.extract_header_year() == 92
  assert pdb_inp.extract_remark_iii_records(iii=2) \
      == ['REMARK   2 RESOLUTION. 1.7  ANGSTROMS.']

def exercise_input_pickling():
  pdb_inp = pdb.input(source_info="file/name", lines=pdb_string_all_sections)
  for p in [pickle, cPickle]:
    s = p.dumps(pdb_inp, 1)
    l = p.loads(s)
    assert not show_diff(l.as_pdb_string(), pdb_inp.as_pdb_string())
    assert l.source_info() == "pickle"
    for section in pdb.input_sections:
      assert not show_diff(
        "\n".join(getattr(l, section)()),
        "\n".join(getattr(pdb_inp, section)()))
    s = "\n".join(l.__getinitargs__()[1])
    d = hashlib_md5(s).hexdigest()
    if (pdb.hierarchy.atom.has_siguij()):
      assert d == "bf987c40cc8672e2f2324d91d6de3e2b"
    else:
      assert d == "7375e96fd52794a785284580730de20c"

def exercise_xray_structure_simple():
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
CRYST1   61.410   54.829   43.543  90.00  90.00  90.00 P 21 21 21    8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.016284  0.000000  0.000000        0.00000
SCALE2      0.000000  0.018239  0.000000        0.00000
SCALE3      0.000000  0.000000  0.022966        0.00000
ATOM      1  N   GLN A   3      35.299  11.075  19.070  1.00 36.89           N
ATOM      2  CA  GLN A   3      34.482   9.927  18.794  0.63 37.88           C
SIGATM    2  CA  GLN A   3       1.200   2.300   3.400  0.04  0.05           C
ANISOU    2  CA  GLN A   3     7794   3221   3376  -1227   1064   2601       C
ATOM      3  Q   GLN A   3      35.130   8.880  17.864  0.84 37.52           C
ANISOU    3  Q   GLN A   3     7875   3041   3340   -981    727   2663       C
SIGUIJ    3  Q   GLN A   3       75     41     40     -1      7     63       C
ATOM      4  O   GLN A   3      34.548   7.819  17.724  1.00 38.54      STUV
ATOM      5 1CB AGLN A   3      32.979  10.223  18.469  1.00 37.80
HETATM    6 CA  AION B   1      32.360  11.092  17.308  0.92 35.96          CA2+
HETATM    7 CA   ION B   2      30.822  10.665  17.190  1.00 36.87
"""))
  for use_scale_matrix_if_available in [False, True]:
    xray_structure = pdb_inp.xray_structure_simple(
      use_scale_matrix_if_available=use_scale_matrix_if_available)
    out = StringIO()
    xray_structure.show_summary(f=out)
    assert not show_diff(out.getvalue(), """\
Number of scatterers: 7
At special positions: 0
Unit cell: (61.41, 54.829, 43.543, 90, 90, 90)
Space group: P 21 21 21 (No. 19)
""")
    out = StringIO()
    xray_structure.show_scatterers(f=out)
    assert not show_diff(out.getvalue(), """\
Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso, Ustar as Uiso
pdb=" N   GLN A   3 " N      4 ( 0.5748  0.2020  0.4380) 1.00 0.4672 [ - ]
pdb=" CA  GLN A   3 " C      4 ( 0.5615  0.1811  0.4316) 0.63 [ - ] 0.4797
     u_cart =  0.779  0.322  0.338 -0.123  0.106  0.260
pdb=" Q   GLN A   3 " C      4 ( 0.5721  0.1620  0.4103) 0.84 [ - ] 0.4752
     u_cart =  0.788  0.304  0.334 -0.098  0.073  0.266
pdb=" O   GLN A   3 " segid="STUV" O      4 ( 0.5626  0.1426  0.4070)\
 1.00 0.4881 [ - ]
pdb="1CB AGLN A   3 " C      4 ( 0.5370  0.1865  0.4242) 1.00 0.4787 [ - ]
pdb="CA  AION B   1 " Ca2+   4 ( 0.5270  0.2023  0.3975) 0.92 0.4554 [ - ]
pdb="CA   ION B   2 " Ca     4 ( 0.5019  0.1945  0.3948) 1.00 0.4670 [ - ]
""")
  #
  xray_structure = pdb_inp.xray_structure_simple(unit_cube_pseudo_crystal=True)
  out = StringIO()
  xray_structure.show_summary(f=out)
  assert not show_diff(out.getvalue(), """\
Number of scatterers: 7
At special positions: 0
Unit cell: (1, 1, 1, 90, 90, 90)
Space group: P 1 (No. 1)
""")
  out = StringIO()
  xray_structure.show_scatterers(f=out)
  assert not show_diff(out.getvalue(), """\
Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso, Ustar as Uiso
pdb=" N   GLN A   3 " N      1 (35.2990 11.0750 19.0700) 1.00 0.4672 [ - ]
pdb=" CA  GLN A   3 " C      1 (34.4820  9.9270 18.7940) 0.63 [ - ] 0.4797
     u_cart =  0.779  0.322  0.338 -0.123  0.106  0.260
pdb=" Q   GLN A   3 " C      1 (35.1300  8.8800 17.8640) 0.84 [ - ] 0.4752
     u_cart =  0.788  0.304  0.334 -0.098  0.073  0.266
pdb=" O   GLN A   3 " segid="STUV" O      1 (34.5480  7.8190 17.7240)\
 1.00 0.4881 [ - ]
pdb="1CB AGLN A   3 " C      1 (32.9790 10.2230 18.4690) 1.00 0.4787 [ - ]
pdb="CA  AION B   1 " Ca2+   1 (32.3600 11.0920 17.3080) 0.92 0.4554 [ - ]
pdb="CA   ION B   2 " Ca     1 (30.8220 10.6650 17.1900) 1.00 0.4670 [ - ]
""")
  #
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM      1  N   GLN A   3      35.299  11.075  99.070  1.00 36.89      STUV A
"""))
  xray_structure = pdb_inp.xray_structure_simple(
    enable_scattering_type_unknown=True)
  assert xray_structure.scatterers()[0].scattering_type == "unknown"
  try:
    pdb_inp.xray_structure_simple()
  except Sorry, e:
    assert not show_diff(str(e), """\
Unknown chemical element type:
  "ATOM      1  N   GLN A   3 .*.STUV A  "
  To resolve this problem, specify a chemical element type in
  columns 77-78 of the PDB file, right justified (e.g. " C").""")
  else: raise Exception_expected
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM      1 1A   GLN A   3      35.299  11.075  99.070  1.00 36.89
"""))
  xray_structure = pdb_inp.xray_structure_simple(
    enable_scattering_type_unknown=True)
  assert xray_structure.scatterers()[0].scattering_type == "unknown"
  try:
    pdb_inp.xray_structure_simple()
  except Sorry, e:
    assert not show_diff(str(e), """\
Unknown chemical element type:
  "ATOM      1 1A   GLN A   3 .*.        "
  To resolve this problem, specify a chemical element type in
  columns 77-78 of the PDB file, right justified (e.g. " C").""")
  else: raise Exception_expected
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM      1  N   GLN A   3      35.299  11.075  99.070  1.00 36.89           Bx5
"""))
  xray_structure = pdb_inp.xray_structure_simple(
    enable_scattering_type_unknown=True)
  assert xray_structure.scatterers()[0].scattering_type == "unknown"
  try:
    pdb_inp.xray_structure_simple()
  except Sorry, e:
    assert not show_diff(str(e), '''\
Unknown charge:
  "ATOM      1  N   GLN A   3 .*.     Bx5"
                                       ^^''')
  else: raise Exception_expected
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM      1  N   GLN A   3      35.299  11.075  99.070  1.00 36.89          Cs3-
"""))
  xray_structure = pdb_inp.xray_structure_simple()
  assert xray_structure.scatterers()[0].scattering_type == "Cs"
  xray_structure = pdb_inp.xray_structure_simple(
    scattering_type_exact=True,
    enable_scattering_type_unknown=True)
  assert xray_structure.scatterers()[0].scattering_type == "unknown"
  out = StringIO()
  xray_structure.scattering_type_registry().show(out=out)
  assert not show_diff(out.getvalue(), """\
Number of scattering types: 1
  Type     Number    sf(0)   Gaussians
   unknown     1      None      None
  sf(0) = scattering factor at diffraction angle 0.
""")
  try:
    pdb_inp.xray_structure_simple(scattering_type_exact=True)
  except Sorry, e:
    assert not show_diff(str(e), '''\
Unknown scattering type:
  "ATOM      1  N   GLN A   3 .*.    Cs3-"
               ^^^^                  ^^^^''')
  else: raise Exception_expected
  #
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
MODEL        1
ATOM      1  N   GLN A   3      35.299  11.075  19.070  1.00 36.89
ATOM      2  CA  GLN A   3      34.482   9.927  18.794  0.63 37.88
ENDMDL
MODEL        2
ATOM      1  N   GLN A   3      25.299   1.075   9.070  0.54 26.89
ATOM      2  CA  GLN A   3      24.482  -1.927   8.794  1.00 27.88
ENDMDL
"""))
  xray_structure = pdb_inp.xray_structure_simple(unit_cube_pseudo_crystal=True)
  out = StringIO()
  xray_structure.show_scatterers(f=out)
  assert not show_diff(out.getvalue(), """\
Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso, Ustar as Uiso
model="   1" pdb=" N   GLN A   3 " N      1 (35.2990 11.0750 19.0700)\
 1.00 0.4672 [ - ]
model="   1" pdb=" CA  GLN A   3 " C      1 (34.4820  9.9270 18.7940)\
 0.63 0.4798 [ - ]
model="   2" pdb=" N   GLN A   3 " N      1 (25.2990  1.0750  9.0700)\
 0.54 0.3406 [ - ]
model="   2" pdb=" CA  GLN A   3 " C      1 (24.4820 -1.9270  8.7940)\
 1.00 0.3531 [ - ]
""")
  xray_structures = pdb_inp.xray_structures_simple(
    unit_cube_pseudo_crystal=True)
  assert len(xray_structures) == 2
  out = StringIO()
  xray_structures[0].show_scatterers(f=out)
  assert not show_diff(out.getvalue(), """\
Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso, Ustar as Uiso
model="   1" pdb=" N   GLN A   3 " N      1 (35.2990 11.0750 19.0700)\
 1.00 0.4672 [ - ]
model="   1" pdb=" CA  GLN A   3 " C      1 (34.4820  9.9270 18.7940)\
 0.63 0.4798 [ - ]
""")
  out = StringIO()
  xray_structures[1].show_scatterers(f=out)
  assert not show_diff(out.getvalue(), """\
Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso, Ustar as Uiso
model="   2" pdb=" N   GLN A   3 " N      1 (25.2990  1.0750  9.0700)\
 0.54 0.3406 [ - ]
model="   2" pdb=" CA  GLN A   3 " C      1 (24.4820 -1.9270  8.7940)\
 1.00 0.3531 [ - ]
""")
  #
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM    369 PEAK PEAK    1      61.114  12.134   8.619  1.00 20.00      PEAK
ATOM    504 SITE SITE    2      67.707   2.505  14.951  1.00 20.00      SITE
"""))
  xray_structure = pdb_inp.xray_structure_simple()
  assert xray_structure.scattering_type_registry().type_index_pairs_as_dict() \
      == {"const": 0}
  assert list(xray_structure.scattering_type_registry().unique_counts) == [2]
  #
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
"""))
  assert pdb_inp.xray_structure_simple().scatterers().size() == 0
  assert len(pdb_inp.xray_structures_simple()) == 1
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
MODEL        1
ENDMDL
"""))
  assert pdb_inp.xray_structure_simple().scatterers().size() == 0
  assert len(pdb_inp.xray_structures_simple()) == 1
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
MODEL        1
ENDMDL
MODEL        2
ENDMDL
"""))
  assert pdb_inp.xray_structure_simple().scatterers().size() == 0
  assert len(pdb_inp.xray_structures_simple()) == 2
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
MODEL        1
ATOM      1  N   GLN A   3      35.299  11.075  99.070  1.00 36.89           O-2
ENDMDL
MODEL        2
ENDMDL
"""))
  assert pdb_inp.xray_structure_simple().scatterers().size() == 1
  xray_structures = pdb_inp.xray_structures_simple()
  assert len(xray_structures) == 2
  assert xray_structures[0].scatterers().size() == 1
  assert xray_structures[1].scatterers().size() == 0
  assert xray_structures[0].scatterers()[0].scattering_type == "O2-"
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
MODEL        1
ENDMDL
MODEL        2
ATOM      1  N   GLN A   3      35.299  11.075  99.070  1.00 36.89          Fe+3
ENDMDL
"""))
  assert pdb_inp.xray_structure_simple().scatterers().size() == 1
  xray_structures = pdb_inp.xray_structures_simple()
  assert len(xray_structures) == 2
  assert xray_structures[0].scatterers().size() == 0
  assert xray_structures[1].scatterers().size() == 1
  assert xray_structures[1].scatterers()[0].scattering_type == "Fe3+"
  input_pdb_string = """\
ATOM      1  N   GLN A   3      35.299  11.075  19.070  1.00 36.89           N
ATOM      2  CA  GLN A   3      34.482   9.927  18.794  0.63 37.88           C0
ATOM      3  C   GLN A   3      35.130   8.880  17.864  0.84 37.52           C 0
ATOM      4  O   GLN A   3      34.548   7.819  17.724  1.00 38.54           O00
ATOM      5 1CB  GLN A   3      32.979  10.223  18.469  1.00 37.80           C 1
HETATM    6 CA   ION B   1      32.360  11.092  17.308  0.92 35.96           X
HETATM    7 CA   ION B   2      30.822  10.665  17.190  1.00 36.87          FE4+
ATOM      8  O   MET A   5       6.215  22.789  24.067  1.00  0.00            -2
"""
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines(input_pdb_string))
  assert [scatterer.scattering_type
    for scatterer in pdb_inp.xray_structure_simple(
      scattering_type_exact=True,
      enable_scattering_type_unknown=True).scatterers()] \
        == ["N", "C", "C", "O", "unknown", "unknown", "unknown", "O2-"]
  #
  cs1 = crystal.symmetry(
    unit_cell=(3.113,3.444,2.572,90,90,90),
    space_group_symbol="P1")
  cs2 = crystal.symmetry(
    unit_cell=(10,20,30,80,85,95),
    space_group_symbol="P-1")
  input_pdb_string = """\
ATOM      1  N   GLN A   3      35.299  11.075  19.070  1.00 36.89           N
ATOM      2  CA  GLN A   3      34.482   9.927  18.794  0.63 37.88           C
"""
  pdb_inp = pdb.input(source_info=None, lines=input_pdb_string)
  xs = pdb_inp.xray_structure_simple()
  assert xs.is_similar_symmetry(cs1)
  xs = pdb_inp.xray_structure_simple(crystal_symmetry=cs2)
  assert xs.is_similar_symmetry(cs2)
  input_pdb_string = """\
CRYST1   10.000   20.000   30.000  80.00  85.00  95.00
ATOM      1  N   GLN A   3      35.299  11.075  19.070  1.00 36.89           N
ATOM      2  CA  GLN A   3      34.482   9.927  18.794  0.63 37.88           C
"""
  pdb_inp = pdb.input(source_info=None, lines=input_pdb_string)
  xs = pdb_inp.xray_structure_simple()
  assert xs.is_similar_symmetry(
    cs1.customized_copy(unit_cell=cs2.unit_cell()))
  xs = pdb_inp.xray_structure_simple(crystal_symmetry=cs2)
  assert xs.is_similar_symmetry(cs2)
  input_pdb_string = """\
CRYST1                                                 P -1
ATOM      1  N   GLN A   3      35.299  11.075  19.070  1.00 36.89           N
ATOM      2  CA  GLN A   3      34.482   9.927  18.794  0.63 37.88           C
"""
  pdb_inp = pdb.input(source_info=None, lines=input_pdb_string)
  xs = pdb_inp.xray_structure_simple()
  assert xs.is_similar_symmetry(cs1.customized_copy(
    space_group_info=cs2.space_group_info()))
  xs = pdb_inp.xray_structure_simple(crystal_symmetry=cs2)
  assert xs.is_similar_symmetry(cs2)
  input_pdb_string = """\
CRYST1   10.000   20.000   30.000  80.00  85.00  95.00 P -1
ATOM      1  N   GLN A   3      35.299  11.075  19.070  1.00 36.89           N
ATOM      2  CA  GLN A   3      34.482   9.927  18.794  0.63 37.88           C
"""
  pdb_inp = pdb.input(source_info=None, lines=input_pdb_string)
  xs = pdb_inp.xray_structure_simple()
  assert xs.is_similar_symmetry(cs2)
  xs = pdb_inp.xray_structure_simple(crystal_symmetry=cs1)
  assert xs.is_similar_symmetry(cs1)

def get_phenix_regression_pdb_file_names():
  pdb_dir = libtbx.env.find_in_repositories("phenix_regression/pdb")
  if (pdb_dir is None): return None
  result = []
  for node in os.listdir(pdb_dir):
    if (not (node.endswith(".pdb") or node.endswith(".ent"))): continue
    result.append(os.path.join(pdb_dir, node))
  assert len(result) != 0
  return result

def exercise(args):
  phenix_regression_pdb_file_names = get_phenix_regression_pdb_file_names()
  forever = "--forever" in args
  while True:
    exercise_hybrid_36()
    exercise_base_256_ordinal()
    exercise_columns_73_76_evaluator(
      pdb_file_names=phenix_regression_pdb_file_names)
    exercise_line_info_exceptions()
    exercise_pdb_input()
    exercise_input_pickling()
    exercise_xray_structure_simple()
    if (not forever): break
  print format_cpu_times()

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
