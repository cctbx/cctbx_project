from iotbx import pdb
import iotbx.pdb.parser
from cctbx.array_family import flex
from scitbx.python_utils import dicts
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal, show_diff
import libtbx.load_env
from cStringIO import StringIO
import sys, os

def exercise_atom():
  a = pdb.atom()
  assert a.name == ""
  a.name = "abcd"
  assert a.name == "abcd"
  try: a.name = "xyzhkl"
  except ValueError, e:
    assert str(e) == "string is too long for name attribute " \
      "(maximum length is 4 characters, 6 given)."
  else: raise RuntimeError("Exception expected.")
  assert a.segid == ""
  a.segid = "stuv"
  assert a.segid == "stuv"
  assert a.element == ""
  a.element = "ca"
  assert a.element == "ca"
  assert a.charge == ""
  a.charge = "2+"
  assert a.charge == "2+"
  assert a.xyz == (0,0,0)
  a.xyz = (1,-2,3)
  assert a.xyz == (1,-2,3)
  assert a.sigxyz == (0,0,0)
  a.sigxyz = (-2,3,1)
  assert a.sigxyz == (-2,3,1)
  assert a.occ == 0
  a.occ = 0.5
  assert a.occ == 0.5
  assert a.sigocc == 0
  a.sigocc = 0.7
  assert a.sigocc == 0.7
  assert a.b == 0
  a.b = 5
  assert a.b == 5
  assert a.sigb == 0
  a.sigb = 7
  assert a.sigb == 7
  assert a.uij == (0,0,0,0,0,0)
  a.uij = (1,-2,3,4,-5,6)
  assert a.uij == (1,-2,3,4,-5,6)
  assert a.siguij == (0,0,0,0,0,0)
  a.siguij = (-2,3,4,-5,6,1)
  assert a.siguij == (-2,3,4,-5,6,1)
  assert not a.hetero
  a.hetero = True
  assert a.hetero
  assert a.tmp == 0
  a.tmp = 3
  assert a.tmp == 3
  #
  a = (pdb.atom()
    .set_name(new_name="NaMe")
    .set_segid(new_segid="sEgI")
    .set_element(new_element="El")
    .set_charge(new_charge="cH")
    .set_xyz(new_xyz=(1.3,2.1,3.2))
    .set_sigxyz(new_sigxyz=(.1,.2,.3))
    .set_occ(new_occ=0.4)
    .set_sigocc(new_sigocc=0.1)
    .set_b(new_b=4.8)
    .set_sigb(new_sigb=0.7)
    .set_uij(new_uij=(1.3,2.1,3.2,4.3,2.7,9.3))
    .set_siguij(new_siguij=(.1,.2,.3,.6,.1,.9))
    .set_hetero(new_hetero=True))
  assert a.name == "NaMe"
  assert a.segid == "sEgI"
  assert a.element == "El"
  assert a.charge == "cH"
  assert approx_equal(a.xyz, (1.3,2.1,3.2))
  assert approx_equal(a.sigxyz, (.1,.2,.3))
  assert approx_equal(a.occ, 0.4)
  assert approx_equal(a.sigocc, 0.1)
  assert approx_equal(a.b, 4.8)
  assert approx_equal(a.sigb, 0.7)
  assert approx_equal(a.uij, (1.3,2.1,3.2,4.3,2.7,9.3))
  assert approx_equal(a.siguij, (.1,.2,.3,.6,.1,.9))
  assert a.hetero
  assert not a.is_alternative()
  assert a.determine_chemical_element_simple() is None
  a.name = "NA  "
  a.element = " N"
  assert a.determine_chemical_element_simple() == " N"
  a.element = "CU"
  assert a.determine_chemical_element_simple() == "CU"
  a.element = "  "
  assert a.determine_chemical_element_simple() == "NA"
  a.name = " D"
  assert a.determine_chemical_element_simple() == " D"
  for d in "0123456789":
    a.name = d+"H"
    assert a.determine_chemical_element_simple() == " H"
  #
  r = pdb.residue()
  ac = pdb.atom(parent=r, other=a)
  assert ac.memory_id() != a.memory_id()
  assert a.parents_size() == 0
  assert ac.parents_size() == 1
  assert ac.parents()[0].memory_id() == r.memory_id()
  assert ac.name == a.name
  assert ac.segid == a.segid
  assert ac.element == a.element
  assert ac.charge == a.charge
  assert ac.xyz == a.xyz
  assert ac.sigxyz == a.sigxyz
  assert ac.occ == a.occ
  assert ac.sigocc == a.sigocc
  assert ac.b == a.b
  assert ac.sigb == a.sigb
  assert ac.uij == a.uij
  assert ac.siguij == a.siguij
  assert ac.hetero == a.hetero
  assert ac.is_alternative() == a.is_alternative()
  assert a.tmp == 0
  #
  r1 = pdb.residue(name="abcd", seq=123, icode="mark")
  r2 = pdb.residue(name="efgh", seq=234, icode="bare")
  assert r1.memory_id() != r2.memory_id()
  a = pdb.atom()
  a.pre_allocate_parents(number_of_additional_parents=2)
  a.add_parent(r1)
  p = a.parents()
  assert len(p) == 1
  assert a.parents_size() == 1
  assert p[0].memory_id() == r1.memory_id()
  a.add_parent(r2)
  p = a.parents()
  assert len(p) == 2
  assert a.parents_size() == 2
  assert p[0].memory_id() == r1.memory_id()
  assert p[1].memory_id() == r2.memory_id()
  del r1
  del p
  p = a.parents()
  assert len(p) == 1
  assert a.parents_size() == 1
  assert p[0].memory_id() == r2.memory_id()
  del r2
  del p
  p = a.parents()
  assert len(p) == 0
  assert a.parents_size() == 0

def exercise_residue():
  r = pdb.residue()
  assert r.name == ""
  assert r.seq == 0
  assert r.icode == ""
  assert r.id() == "       0"
  assert r.link_to_previous
  r = pdb.residue(name="xyzw", seq=123, icode="ijkl", link_to_previous=False)
  assert r.name == "xyzw"
  assert r.seq == 123
  assert r.icode == "ijkl"
  assert r.id() == "xyzw 123ijkl"
  assert not r.link_to_previous
  r.link_to_previous = True
  assert r.link_to_previous
  r.name = "foo"
  r.seq = -3
  r.icode = "bar"
  assert r.id() == "foo   -3bar"
  #
  f = pdb.conformer(id="a")
  r.add_atom(new_atom=pdb.atom().set_name(new_name="n"))
  rc = pdb.residue(parent=f, other=r)
  assert rc.memory_id() != r.memory_id()
  assert r.parent() is None
  assert rc.parent().memory_id() == f.memory_id()
  assert rc.id() == r.id()
  assert rc.link_to_previous
  assert rc.atoms_size() == 1
  assert rc.atoms()[0].memory_id() != r.atoms()[0].memory_id()
  assert rc.atoms()[0].name == "n"
  r.add_atom(new_atom=pdb.atom().set_name(new_name="o"))
  assert rc.atoms_size() == 1
  assert r.atoms_size() == 2
  #
  c1 = pdb.conformer(id="a")
  c2 = pdb.conformer(id="b")
  assert c1.memory_id() != c2.memory_id()
  r = pdb.residue(parent=c1)
  assert r.parent().memory_id() == c1.memory_id()
  r = pdb.residue()
  assert r.parent() is None
  r.set_parent(new_parent=c1)
  assert r.parent().memory_id() == c1.memory_id()
  r.set_parent(new_parent=c2)
  assert r.parent().memory_id() == c2.memory_id()
  del c2
  assert r.parent() is None
  assert not r.parent_chain_has_multiple_conformers()
  #
  c1 = pdb.conformer(id="a")
  r = c1.new_residue(name="a", seq=1, icode="i", link_to_previous=False)
  assert r.parent().memory_id() == c1.memory_id()
  del c1
  assert r.parent() is None
  #
  r.pre_allocate_atoms(number_of_additional_atoms=2)
  assert r.atoms_size() == 0
  assert len(r.atoms()) == 0
  r.add_atom(new_atom=pdb.atom().set_name(new_name="ca"))
  assert r.atoms_size() == 1
  assert len(r.atoms()) == 1
  r.add_atom(new_atom=pdb.atom().set_name(new_name="n"))
  assert r.atoms_size() == 2
  assert len(r.atoms()) == 2
  assert [atom.name for atom in r.atoms()] == ["ca", "n"]
  r.new_atoms(number_of_additional_atoms=3)
  assert r.atoms_size() == 5
  assert len(r.atoms()) == 5
  for atom in r.atoms():
    assert atom.parents()[0].memory_id() == r.memory_id()
  assert r.number_of_alternative_atoms() == 0
  for atom in r.atoms():
    assert atom.tmp == 0
  assert r.reset_atom_tmp(new_value=7) == 5
  for atom in r.atoms():
    assert atom.tmp == 7
  assert r.center_of_geometry() == (0,0,0)

def exercise_conformer():
  c = pdb.conformer()
  assert c.id == ""
  c = pdb.conformer(id="a")
  assert c.id == "a"
  c.id = "x"
  assert c.id == "x"
  #
  c1 = pdb.chain(id="a")
  c2 = pdb.chain(id="b")
  assert c1.memory_id() != c2.memory_id()
  f = pdb.conformer(parent=c1)
  assert f.parent().memory_id() == c1.memory_id()
  f = pdb.conformer()
  assert f.parent() is None
  f.set_parent(new_parent=c1)
  assert f.parent().memory_id() == c1.memory_id()
  f.set_parent(new_parent=c2)
  assert f.parent().memory_id() == c2.memory_id()
  del c2
  assert f.parent() is None
  assert not f.parent_chain_has_multiple_conformers()
  #
  c1 = pdb.conformer(id="a")
  c1.pre_allocate_residues(number_of_additional_residues=2)
  assert c1.residues_size() == 0
  assert len(c1.residues()) == 0
  c1.new_residues(number_of_additional_residues=2)
  assert c1.residues_size() == 2
  assert len(c1.residues()) == 2
  for residue in c1.residues():
    assert residue.parent().memory_id() == c1.memory_id()
  assert c1.number_of_atoms() == 0
  assert c1.number_of_alternative_atoms() == 0
  assert c1.reset_atom_tmp(new_value=8) == 0
  assert list(c1.residue_centers_of_geometry()) == [(0,0,0), (0,0,0)]

def exercise_chain():
  c = pdb.chain()
  assert c.id == ""
  c = pdb.chain(id="a")
  assert c.id == "a"
  c.id = "x"
  assert c.id == "x"
  #
  m1 = pdb.model(id=1)
  m2 = pdb.model(id=2)
  assert m1.memory_id() != m2.memory_id()
  c = pdb.chain(parent=m1)
  assert c.parent().memory_id() == m1.memory_id()
  c = pdb.chain()
  assert c.parent() is None
  c.set_parent(new_parent=m1)
  assert c.parent().memory_id() == m1.memory_id()
  c.set_parent(new_parent=m2)
  assert c.parent().memory_id() == m2.memory_id()
  del m2
  assert c.parent() is None
  assert not c.has_multiple_conformers()
  #
  c = pdb.chain()
  c.pre_allocate_conformers(number_of_additional_conformers=2)
  assert c.conformers_size() == 0
  assert len(c.conformers()) == 0
  c.new_conformers(number_of_additional_conformers=2)
  assert c.conformers_size() == 2
  assert len(c.conformers()) == 2
  for conformer in c.conformers():
    assert conformer.parent().memory_id() == c.memory_id()
  assert c.reset_atom_tmp(new_value=9) == 0
  assert c.extract_sites_cart().size() == 0

def exercise_model():
  m = pdb.model()
  assert m.id == 0
  m = pdb.model(id=42)
  assert m.id == 42
  m.id = -23
  assert m.id == -23
  #
  m = pdb.model(id=1)
  assert m.parent() is None
  m.pre_allocate_chains(number_of_additional_chains=2)
  assert m.chains_size() == 0
  assert len(m.chains()) == 0
  ch_a = m.new_chain(chain_id="a")
  assert ch_a.parent().memory_id() == m.memory_id()
  assert m.chains_size() == 1
  assert len(m.chains()) == 1
  ch_b = pdb.chain(id="b")
  assert ch_b.parent() is None
  m.adopt_chain(new_chain=ch_b)
  assert m.chains_size() == 2
  chains = m.chains()
  assert len(chains) == 2
  assert chains[0].memory_id() == ch_a.memory_id()
  assert chains[1].memory_id() == ch_b.memory_id()
  m.new_chains(number_of_additional_chains=3)
  assert m.chains_size() == 5
  assert len(m.chains()) == 5
  for chain in m.chains():
    assert chain.parent().memory_id() == m.memory_id()
  assert m.reset_atom_tmp(new_value=2) == 0

def exercise_hierarchy():
  h = pdb.hierarchy()
  m = pdb.model()
  m.set_parent(new_parent=h)
  assert m.parent().memory_id() == h.memory_id()
  m = pdb.model(parent=h)
  assert m.parent().memory_id() == h.memory_id()
  assert m.id == 0
  m = pdb.model(parent=h, id=2)
  assert m.parent().memory_id() == h.memory_id()
  assert m.id == 2
  del h
  assert m.parent() is None
  #
  h = pdb.hierarchy()
  assert h.info.size() == 0
  h.info.append("abc")
  assert h.info.size() == 1
  h.info = flex.std_string(["a", "b"])
  assert h.info.size() == 2
  h.pre_allocate_models(number_of_additional_models=2)
  assert h.models_size() == 0
  assert len(h.models()) == 0
  m_a = h.new_model(model_id=3)
  assert m_a.parent().memory_id() == h.memory_id()
  assert h.models_size() == 1
  assert len(h.models()) == 1
  m_b = pdb.model(id=5)
  assert m_b.parent() is None
  h.adopt_model(new_model=m_b)
  assert h.models_size() == 2
  models = h.models()
  assert len(models) == 2
  assert models[0].memory_id() == m_a.memory_id()
  assert models[1].memory_id() == m_b.memory_id()
  h.new_models(number_of_additional_models=3)
  assert h.models_size() == 5
  assert len(h.models()) == 5
  for model in h.models():
    assert model.parent().memory_id() == h.memory_id()
  h.reset_atom_tmp(new_value=1) == 0

def check_hierarchy(
      hierarchy,
      expected_formatted=None,
      expected_overall_counts=None):
  out = StringIO()
  hierarchy.show(out=out)
  if (expected_formatted is None or expected_formatted == "None\n"):
    sys.stdout.write(out.getvalue())
    print "#"*79
  else:
    assert not show_diff(out.getvalue(), expected_formatted)
  overall_counts = hierarchy.overall_counts()
  if (expected_overall_counts is not None):
    assert overall_counts.chain_ids == expected_overall_counts.chain_ids
    assert overall_counts.conformer_ids ==expected_overall_counts.conformer_ids
    assert overall_counts.residue_names ==expected_overall_counts.residue_names
    assert overall_counts.residue_name_classes \
        == expected_overall_counts.residue_name_classes
  assert overall_counts.n_chains \
      == sum(overall_counts.chain_ids.values())
  assert overall_counts.n_conformers \
      == sum(overall_counts.conformer_ids.values())
  if (overall_counts.n_conformers == overall_counts.n_chains):
    assert overall_counts.n_residues \
        == sum(overall_counts.residue_names.values())
  assert sum(overall_counts.residue_name_classes.values()) \
    == sum(overall_counts.residue_names.values())

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
  #
  lines = flex.split_lines("""\
HEADER    HYDROLASE(METALLOPROTEINASE)            17-NOV-93   1THL
ATOM      1  N   ILE     1       9.581  51.813  -0.720  1.00 31.90      1THL 158
ATOM      2  CA  ILE     1       8.335  52.235  -0.041  1.00 52.95      1THL 159
ATOM      3  C   ILE     1       7.959  53.741   0.036  1.00 26.88      1THL 160
END
""")
  cpp_eval = pdb.columns_73_76_evaluator(lines=lines)
  assert cpp_eval.finding == "Exactly one common label in columns 73-76."
  assert cpp_eval.is_old_style

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
    assert pdb_inp.input_atom_labels_list().size() == 0
    assert pdb_inp.atom_serial_number_strings().size() == 0
    assert pdb_inp.atoms().size() == 0
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
    out = StringIO()
    pdb_inp.construct_hierarchy().show(out=out)
    assert len(out.getvalue()) == 0
    assert pdb_inp.number_of_chains_with_altloc_mix() == 0
    assert pdb_inp.i_seqs_alternative_group_with_blank_altloc().size() == 0
    assert pdb_inp.i_seqs_alternative_group_without_blank_altloc().size() == 0
    assert pdb_inp.atom_element_counts() == {}
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
HETATM    3  C   MET A   2       7.478  21.387  22.491  1.00  0.00           C
ATOM      4  O   MET A   2       8.406  20.895  23.132  1.00  0.00           O
ENDMDL
MODEL 3
HETATM    9 2H3  MPR B   5      16.388   0.289   6.613  1.00  0.08
SIGATM    9 2H3  MPR B   5       0.155   0.175   0.155  0.00  0.05
ANISOU    9 2H3  MPR B   5      848    848    848      0      0      0
SIGUIJ    9 2H3  MPR B   5      510    510    510      0      0      0
TER
ATOM     10  N   CYS C   6      14.270   2.464   3.364  1.00  0.07
SIGATM   10  N   CYS C   6       0.012   0.012   0.011  0.00  0.00
ANISOU   10  N   CYS C   6      788    626    677   -344    621   -232
SIGUIJ   10  N   CYS C   6        3     13      4     11      6     13
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
    assert list(pdb_inp.atom_serial_number_strings()) \
        == ["    1", "    2", "    3", "    4", "    9", "   10"]
    assert [atom.element for atom in pdb_inp.atoms()] \
        == [" N", " C", " C", " O", "  ", "  "]
    assert [atom.parents() for atom in pdb_inp.atoms()] \
        == [[], [], [], [], [], []]
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
    assert pdb_inp.chain_selection_cache().keys() == ["A", "B", "C"]
    assert [list(v) for v in pdb_inp.chain_selection_cache().values()] \
        == [[0,1,2,3], [4], [5]]
    for resseq,i_seqs in [(1,[0,1]),(2,[2,3]),(5,[4]),(6,[5])]:
      assert list(pdb_inp.resseq_selection_cache()[resseq+999]) == i_seqs
    for i,i_seqs in enumerate(pdb_inp.resseq_selection_cache()):
      if (i not in [j+999 for j in [1,2,5,6]]):
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
    check_hierarchy(
      hierarchy=pdb_inp.construct_hierarchy(),
      expected_formatted="""\
model id=1 #chains=1
  chain id="A" #conformers=1
    conformer id=" " #residues=2 #atoms=4
      residue name="MET " seq=   1 icode=" " #atoms=2
         " N  "
         " CA "
      ### chain break ###
      residue name="MET " seq=   2 icode=" " #atoms=2
         " C  "
         " O  "
model id=3 #chains=2
  chain id="B" #conformers=1
    conformer id=" " #residues=1 #atoms=1
      residue name="MPR " seq=   5 icode=" " #atoms=1
         "2H3 "
  chain id="C" #conformers=1
    conformer id=" " #residues=1 #atoms=1
      residue name="CYS " seq=   6 icode=" " #atoms=1
         " N  "
""",
      expected_overall_counts=dicts.easy(
        chain_ids={"A": 1, "B": 1, "C": 1},
        conformer_ids={" ": 3},
        residue_names={"MPR ": 1, "MET ": 2, "CYS ": 1},
        residue_name_classes={"other": 1, "common_amino_acid": 3}))
  #
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM      1  N   MET A   1       6.215  22.789  24.067  1.00  0.00           N
ATOM      2  N   MET A   1       2.615  27.289  20.467  1.00  0.00           O
"""))
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=1
    conformer id=" " #residues=1 #atoms=2
      residue name="MET " seq=   1 icode=" " #atoms=2
         " N  "
         " N  "
""")
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
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=1 #chains=1
  chain id="A" #conformers=1
    conformer id=" " #residues=1 #atoms=2
      residue name="MET " seq=   1 icode=" " #atoms=2
         " N  "
         " N  "
model id=2 #chains=1
  chain id="A" #conformers=1
    conformer id=" " #residues=1 #atoms=2
      residue name="MET " seq=   1 icode=" " #atoms=2
         " N  "
         " N  "
""")
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
BREAK
ATOM      6  C   MET A   2       2.615  27.289  20.467  1.00  0.00           O
ATOM      7  C   MET A   1       2.615  27.289  20.467  1.00  0.00           O
ATOM      8  C   MET A   1       2.615  27.289  20.467  1.00  0.00           O
"""))
  hierarchy = pdb_inp.construct_hierarchy()
  check_hierarchy(
    hierarchy=hierarchy,
    expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=1
    conformer id=" " #residues=5 #atoms=8
      residue name="MET " seq=   1 icode=" " #atoms=1
         " N  "
      residue name="MET " seq=   2 icode=" " #atoms=1
         " N  "
      residue name="MET " seq=   1 icode=" " #atoms=3
         " N  "
         " N  "
         " C  "
      ### chain break ###
      residue name="MET " seq=   2 icode=" " #atoms=1
         " C  "
      residue name="MET " seq=   1 icode=" " #atoms=2
         " C  "
         " C  "
""")
  c = hierarchy.models()[0].chains()[0]
  sites_cart = c.extract_sites_cart()
  assert sites_cart.size() == 8
  assert approx_equal(sites_cart.mean(), [3.515, 26.164, 21.3669999])
  f = c.conformers()[0]
  for r,g in zip(f.residues(), f.residue_centers_of_geometry()):
    assert approx_equal(r.center_of_geometry(), g)
  assert approx_equal(f.residue_centers_of_geometry()[2],(3.815,25.789,21.667))
  for i,r in enumerate(f.residues()):
    if (i in [0, 3]):
      assert not r.link_to_previous
    else:
      assert r.link_to_previous
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
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=2
  chain id=" " #conformers=1
    conformer id=" " #residues=1 #atoms=2
      residue name="LYS " seq= 109 icode=" " #atoms=2
         " CB "
         " CG "
  chain id=" " #conformers=1
    conformer id=" " #residues=1 #atoms=2
      residue name="LYS " seq= 109 icode=" " #atoms=2
         " CB "
         " CG "
""")
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
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id="C" #conformers=1
    conformer id="A" #residues=1 #atoms=2
      residue name="R123" seq=1234 icode="I" #atoms=2
         "N123"
         "N123"
""")
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
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id="C" #conformers=1
    conformer id="A" #residues=1 #atoms=2
      residue name="R123" seq=1234 icode="I" #atoms=2
         "N123"
         "N123"
""")
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
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id=" " #conformers=1
    conformer id=" " #residues=1 #atoms=2
      residue name="    " seq=   0 icode=" " #atoms=2
         "    "
         "    "
""")
  #
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
"""))
  assert list(pdb_inp.model_indices()) == []
  assert list(pdb_inp.chain_indices()) == []
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
""",
    expected_overall_counts=dicts.easy(
      chain_ids={},
      conformer_ids={},
      residue_names={},
      residue_name_classes={}))
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM
"""))
  assert list(pdb_inp.model_indices()) == [1]
  assert [list(v) for v in pdb_inp.chain_indices()] == [[1]]
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id=" " #conformers=1
    conformer id=" " #residues=1 #atoms=1
      residue name="    " seq=   0 icode=" " #atoms=1
         "    "
""")
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
MODEL        1
ENDMDL
"""))
  assert list(pdb_inp.model_indices()) == [0]
  assert [list(v) for v in pdb_inp.chain_indices()] == [[]]
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=1 #chains=0
""")
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
MODEL        1
ATOM
ENDMDL
"""))
  assert list(pdb_inp.model_indices()) == [1]
  assert [list(v) for v in pdb_inp.chain_indices()] == [[1]]
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=1 #chains=1
  chain id=" " #conformers=1
    conformer id=" " #residues=1 #atoms=1
      residue name="    " seq=   0 icode=" " #atoms=1
         "    "
""")
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
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=1 #chains=0
model id=2 #chains=0
""")
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
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=1 #chains=0
model id=2 #chains=1
  chain id=" " #conformers=1
    conformer id=" " #residues=1 #atoms=1
      residue name="    " seq=   0 icode=" " #atoms=1
         "    "
""")
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
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=1 #chains=3
  chain id="C" #conformers=1
    conformer id=" " #residues=1 #atoms=2
      residue name="    " seq=   0 icode=" " #atoms=2
         "    "
         "    "
  chain id="D" #conformers=1
    conformer id=" " #residues=1 #atoms=1
      residue name="    " seq=   0 icode=" " #atoms=1
         "    "
  chain id="E" #conformers=1
    conformer id=" " #residues=1 #atoms=2
      residue name="    " seq=   0 icode=" " #atoms=2
         "    "
         "    "
model id=2 #chains=3
  chain id="C" #conformers=1
    conformer id=" " #residues=1 #atoms=2
      residue name="    " seq=   0 icode=" " #atoms=2
         "    "
         "    "
  chain id="D" #conformers=1
    conformer id=" " #residues=1 #atoms=2
      residue name="    " seq=   0 icode=" " #atoms=2
         "    "
         "    "
  chain id="E" #conformers=1
    conformer id=" " #residues=1 #atoms=1
      residue name="    " seq=   0 icode=" " #atoms=1
         "    "
model id=3 #chains=3
  chain id=" " #conformers=1
    conformer id=" " #residues=1 #atoms=1
      residue name="    " seq=   0 icode=" " #atoms=1
         "    "
  chain id=" " #conformers=1
    conformer id=" " #residues=1 #atoms=2
      residue name="    " seq=   0 icode=" " #atoms=2
         "    "
         "    "
  chain id="E" #conformers=1
    conformer id=" " #residues=1 #atoms=2
      residue name="    " seq=   0 icode=" " #atoms=2
         "    "
         "    "
model id=4 #chains=2
  chain id=" " #conformers=1
    conformer id=" " #residues=1 #atoms=3
      residue name="    " seq=   0 icode=" " #atoms=3
         "    "
         "    "
         "    "
  chain id="E" #conformers=1
    conformer id=" " #residues=1 #atoms=2
      residue name="    " seq=   0 icode=" " #atoms=2
         "    "
         "    "
""",
    expected_overall_counts=dicts.easy(
      chain_ids={"C": 2, "D": 2, "E": 4, " ": 3},
      conformer_ids={" ": 11},
      residue_names={"    ": 11},
      residue_name_classes={"other": 11}))
  #
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM     54  CA  GLY A   9
"""))
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=1
    conformer id=" " #residues=1 #atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=1
         " CA "
""")
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM     54  CA BGLY A   9
"""))
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=1
    conformer id="B" #residues=1 #atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=1
         " CA "
""")
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM     54  CA BGLY A   9
ATOM     55  CA CGLY A   9
"""))
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=2
    conformer id="B" #residues=1 #atoms=1 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=1
       & " CA "
    conformer id="C" #residues=1 #atoms=1 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=1
       & " CA "
""")
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM     54  CA  GLY A   9
ATOM     55  CA CGLY A   9
"""))
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=2
    conformer id=" " #residues=1 #atoms=1 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=1
       & " CA "
    conformer id="C" #residues=1 #atoms=1 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=1
       & " CA "
""")
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM     54  CA BGLY A   9
ATOM     55  CA  GLY A   9
"""))
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=2
    conformer id=" " #residues=1 #atoms=1 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=1
       & " CA "
    conformer id="B" #residues=1 #atoms=1 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=1
       & " CA "
""")
  perm = flex.size_t((0,1,2))
  for i_trial in xrange(6):
    pdb_inp = pdb.input(
      source_info=None,
      lines=flex.split_lines("""\
ATOM     54  CA BGLY A   9
ATOM     55  CA  GLY A   9
ATOM     56  CA AGLY A   9
""").select(perm))
    check_hierarchy(
      hierarchy=pdb_inp.construct_hierarchy(),
      expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=3
    conformer id=" " #residues=1 #atoms=1 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=1
       & " CA "
    conformer id="A" #residues=1 #atoms=1 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=1
       & " CA "
    conformer id="B" #residues=1 #atoms=1 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=1
       & " CA "
""")
    perm.next_permutation()
  perm = flex.size_t((0,1,2))
  for i_trial in xrange(6):
    pdb_inp = pdb.input(
      source_info=None,
      lines=flex.split_lines("""\
ATOM     54  CA BGLY A   9
ATOM     55  CA  GLY A   9
ATOM     56  CA AGLY A   9
ATOM     57  O   GLY A   9
""").select(perm.concatenate(flex.size_t([3]))))
    check_hierarchy(
      hierarchy=pdb_inp.construct_hierarchy(),
      expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=3
    conformer id=" " #residues=1 #atoms=2 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=2
       & " CA "
         " O  "
    conformer id="A" #residues=1 #atoms=2 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=2
       & " CA "
         " O  "
    conformer id="B" #residues=1 #atoms=2 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=2
       & " CA "
         " O  "
""")
  perm = flex.size_t((1,2,3))
  for i_trial in xrange(6):
    pdb_inp = pdb.input(
      source_info=None,
      lines=flex.split_lines("""\
ATOM     53  O   GLY A   9
ATOM     54  CA BGLY A   9
ATOM     55  CA  GLY A   9
ATOM     56  CA AGLY A   9
""").select(flex.size_t([0]).concatenate(perm)))
    check_hierarchy(
      hierarchy=pdb_inp.construct_hierarchy(),
      expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=3
    conformer id=" " #residues=1 #atoms=2 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=2
         " O  "
       & " CA "
    conformer id="A" #residues=1 #atoms=2 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=2
         " O  "
       & " CA "
    conformer id="B" #residues=1 #atoms=2 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=2
         " O  "
       & " CA "
""")
    perm.next_permutation()
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM     53  CA BGLY A   9
ATOM     54  O   GLY A   9
ATOM     55  CA  GLY A   9
ATOM     56  CA AGLY A   9
"""))
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=3
    conformer id=" " #residues=1 #atoms=2 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=2
         " O  "
       & " CA "
    conformer id="A" #residues=1 #atoms=2 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=2
         " O  "
       & " CA "
    conformer id="B" #residues=1 #atoms=2 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=2
       & " CA "
         " O  "
""")
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM     54  CA  GLY A   9
ATOM     55  CA  GLY A   9
ATOM     56  CA BGLY A   9
"""))
  assert [list(a) for a in pdb_inp.find_duplicate_atom_labels()] == [[0, 1]]
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=2
    conformer id=" " #residues=1 #atoms=2 #alt. atoms=2
      residue name="GLY " seq=   9 icode=" " #atoms=2
       & " CA "
       & " CA "
    conformer id="B" #residues=1 #atoms=1 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=1
       & " CA "
""")
  assert pdb_inp.number_of_chains_with_altloc_mix() == 0
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM     54  CA  GLY A   9
ATOM     55  CA BGLY A   9
ATOM     56  CA  GLY A   9
"""))
  assert [list(a) for a in pdb_inp.find_duplicate_atom_labels()] == [[0, 2]]
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=2
    conformer id=" " #residues=1 #atoms=2 #alt. atoms=2
      residue name="GLY " seq=   9 icode=" " #atoms=2
       & " CA "
       & " CA "
    conformer id="B" #residues=1 #atoms=1 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=1
       & " CA "
""")
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM     54  CA BGLY A   9
ATOM     55  CA  GLY A   9
ATOM     56  CA  GLY A   9
"""))
  assert [list(a) for a in pdb_inp.find_duplicate_atom_labels()] == [[1, 2]]
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=2
    conformer id=" " #residues=1 #atoms=2 #alt. atoms=2
      residue name="GLY " seq=   9 icode=" " #atoms=2
       & " CA "
       & " CA "
    conformer id="B" #residues=1 #atoms=1 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=1
       & " CA "
""")
  perm = flex.size_t((0,1,2))
  for i_trial in xrange(6):
    pdb_inp = pdb.input(
      source_info=None,
      lines=flex.split_lines("""\
ATOM     54  CA  GLY A   9
ATOM     55  CA BGLY A   9
ATOM     56  CA BGLY A   9
""").select(perm))
    invp = perm.inverse_permutation()
    assert [list(a) for a in pdb_inp.find_duplicate_atom_labels()] \
        == [sorted([invp[1], invp[2]])]
    check_hierarchy(
      hierarchy=pdb_inp.construct_hierarchy(),
      expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=2
    conformer id=" " #residues=1 #atoms=1 #alt. atoms=1
      residue name="GLY " seq=   9 icode=" " #atoms=1
       & " CA "
    conformer id="B" #residues=1 #atoms=2 #alt. atoms=2
      residue name="GLY " seq=   9 icode=" " #atoms=2
       & " CA "
       & " CA "
""")
    perm.next_permutation()
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
HEADER    SUGAR BINDING PROTEIN                   20-MAR-00   1EN2
ATOM     54  CA  GLY A   9
ATOM     58  CA AGLY A  10
ATOM     62  CA BSER A  10
ATOM     68  CA  THR A  11
ATOM     75  CA  CYS A  12
ATOM     81  CA  PRO A  13
ATOM     88  CA AALA A  14
ATOM     93  CA BGLY A  14
ATOM     97  CA  LEU A  15
ATOM    108  CA ATRP A  16
ATOM    122  CA BARG A  16
ATOM    140  CA  CYS A  17
"""))
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id="A" #conformers=2
    conformer id="A" #residues=9 #atoms=9 #alt. atoms=3
      residue name="GLY " seq=   9 icode=" " #atoms=1
         " CA "
      residue name="GLY " seq=  10 icode=" " #atoms=1
       & " CA "
      residue name="THR " seq=  11 icode=" " #atoms=1
         " CA "
      residue name="CYS " seq=  12 icode=" " #atoms=1
         " CA "
      residue name="PRO " seq=  13 icode=" " #atoms=1
         " CA "
      residue name="ALA " seq=  14 icode=" " #atoms=1
       & " CA "
      residue name="LEU " seq=  15 icode=" " #atoms=1
         " CA "
      residue name="TRP " seq=  16 icode=" " #atoms=1
       & " CA "
      residue name="CYS " seq=  17 icode=" " #atoms=1
         " CA "
    conformer id="B" #residues=9 #atoms=9 #alt. atoms=3
      residue name="GLY " seq=   9 icode=" " #atoms=1
         " CA "
      residue name="SER " seq=  10 icode=" " #atoms=1
       & " CA "
      residue name="THR " seq=  11 icode=" " #atoms=1
         " CA "
      residue name="CYS " seq=  12 icode=" " #atoms=1
         " CA "
      residue name="PRO " seq=  13 icode=" " #atoms=1
         " CA "
      residue name="GLY " seq=  14 icode=" " #atoms=1
       & " CA "
      residue name="LEU " seq=  15 icode=" " #atoms=1
         " CA "
      residue name="ARG " seq=  16 icode=" " #atoms=1
       & " CA "
      residue name="CYS " seq=  17 icode=" " #atoms=1
         " CA "
""",
    expected_overall_counts=dicts.easy(
      chain_ids={"A": 1},
      conformer_ids={"A": 1, "B": 1},
      residue_names={
        "ALA ": 1, "LEU ": 1, "TRP ": 1, "CYS ": 2, "THR ": 1, "GLY ": 3,
        "SER ": 1, "ARG ": 1, "PRO ": 1},
      residue_name_classes={"common_amino_acid": 12}))
  assert pdb_inp.number_of_chains_with_altloc_mix() == 0
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
HEADER    HYDROLASE(ENDORIBONUCLEASE)             07-JUN-93   1RGA
ATOM    247  OG  SER    35
ATOM    242  N   SER    35
ATOM    244  C   SER    35
ATOM    243  CA  SER    35
ATOM    248  OG BSER    35
ATOM    246  CB BSER    35
ATOM    245  O   SER    35
ATOM    246  CB  SER    35
BREAK
ATOM    250  N  BASN    36
ATOM    252  CA BASN    36
ATOM    253  C   ASN    36
ATOM    251  CA  ASN    36
ATOM    249  N   ASN    36
ATOM    256  O  BASN    36
ATOM    257  CB  ASN    36
ATOM    255  O   ASN    36
ATOM    258  CB BASN    36
"""))
  hierarchy = pdb_inp.construct_hierarchy()
  check_hierarchy(
    hierarchy=hierarchy,
    expected_formatted="""\
model id=0 #chains=1
  chain id=" " #conformers=2
    conformer id=" " #residues=2 #atoms=11 #alt. atoms=6
      residue name="SER " seq=  35 icode=" " #atoms=6
       & " OG "
         " N  "
         " C  "
         " CA "
         " O  "
       & " CB "
      ### chain break ###
      residue name="ASN " seq=  36 icode=" " #atoms=5
         " C  "
       & " CA "
       & " N  "
       & " CB "
       & " O  "
    conformer id="B" #residues=2 #atoms=11 #alt. atoms=6
      residue name="SER " seq=  35 icode=" " #atoms=6
         " N  "
         " C  "
         " CA "
       & " OG "
       & " CB "
         " O  "
      ### chain break ###
      residue name="ASN " seq=  36 icode=" " #atoms=5
       & " N  "
       & " CA "
         " C  "
       & " O  "
       & " CB "
""")
  assert pdb_inp.number_of_chains_with_altloc_mix() == 0
  c = hierarchy.models()[0].chains()[0]
  for f in c.conformers():
    assert [r.number_of_alternative_atoms() for r in f.residues()] == [2,4]
    assert f.number_of_atoms() == 11
    assert f.number_of_alternative_atoms() == 6
  c.reset_atom_tmp(new_value=1)
  assert c.reset_atom_tmp(new_value=0) == 17
  try: pdb_inp.construct_hierarchy()
  except RuntimeError, e:
    assert str(e).endswith(
      "): SCITBX_ASSERT(!construct_hierarchy_was_called_before) failure.")
  else: raise RuntimeError("Exception expected.")
  #
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
HEADER    OXIDOREDUCTASE                          25-NOV-00   1HEU
ATOM      2  CA ASER A   1
ATOM      5  CB ASER A   1
ATOM      8  CA BSER A   1
ATOM     11  CB BSER A   1
ATOM   2092  CA  GLU A 256
ATOM   2095  CB  GLU A 256
ATOM   2100  CB AGLU A 256
ATOM   2105  CB BGLU A 256
END
"""))
  pdb_inp.construct_hierarchy()
  assert pdb_inp.number_of_chains_with_altloc_mix() == 1
  assert list(pdb_inp.i_seqs_alternative_group_with_blank_altloc()) == [5,6,7]
  assert list(pdb_inp.i_seqs_alternative_group_without_blank_altloc()) == [0,2]
  assert pdb_inp.atom_element_counts() == {"  ": 8}
  assert not show_diff(pdb_inp.have_altloc_mix_message(prefix="="), '''\
=mix of alternative groups with and without blank altlocs:
=  alternative group with blank altloc:
=    " CB  GLU A 256 "
=    " CB AGLU A 256 "
=    " CB BGLU A 256 "
=  alternative group without blank altloc:
=    " CA ASER A   1 "
=    " CA BSER A   1 "''')
  try: pdb_inp.raise_altloc_mix_if_necessary()
  except RuntimeError, e:
    assert not show_diff(str(e), '''\
mix of alternative groups with and without blank altlocs:
  alternative group with blank altloc:
    " CB  GLU A 256 "
    " CB AGLU A 256 "
    " CB BGLU A 256 "
  alternative group without blank altloc:
    " CA ASER A   1 "
    " CA BSER A   1 "''')
  else: raise RuntimeError("Exception expected.")
  #
  try: pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
REMARK
ATOM      1  CB  LYS   109
BREAK
ATOM      2  CG  LYS   109
""")).construct_hierarchy()
  except RuntimeError, e:
    assert not show_diff(str(e), "Misplaced BREAK record (input line 3).")
  else: raise RuntimeError("Exception expected.")
  try: pdb.input(
    source_info="file abc",
    lines=flex.split_lines("""\
REMARK
ATOM      1  CA  LYS   109
ATOM      2  CB  LYS   109
BREAK
ATOM      3  CA  LYS   110
BREAK
ATOM      4  CB  LYS   110
""")).construct_hierarchy()
  except RuntimeError, e:
    assert not show_diff(str(e), "Misplaced BREAK record (file abc, line 6).")
  else: raise RuntimeError("Exception expected.")
  #
  open("tmp.pdb", "w")
  pdb_inp = pdb.input(file_name="tmp.pdb")
  assert pdb_inp.source_info() == "file tmp.pdb"
  open("tmp.pdb", "w").write("""\
ATOM      1  CA  SER     1       1.212 -12.134   3.757  1.00  0.00
ATOM      2  CA  LEU     2       1.118  -9.777   0.735  1.00  0.00
""")
  pdb_inp = pdb.input(file_name="tmp.pdb")
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id=" " #conformers=1
    conformer id=" " #residues=2 #atoms=2
      residue name="SER " seq=   1 icode=" " #atoms=1
         " CA "
      residue name="LEU " seq=   2 icode=" " #atoms=1
         " CA "
""")
  try: pdb.input(file_name="")
  except RuntimeError, e:
    assert not show_diff(str(e), 'Cannot open file for reading: ""')
  else: raise RuntimeError("Exception expected.")
  #
  out = StringIO()
  pdb_inp, hierarchy = pdb.show_summary(
    source_info="s",
    lines=flex.split_lines("""\
MODEL 1
ATOM      1  N   MET A   1
ATOM      2  N   MET A   1
ENDMDL
MODEL 2
ATOM      1  N  AMET A   1
ATOM      2  N  BMET A   1
ENDMDL
"""),
    level_id="atom",
    out=out,
    prefix="%#")
  assert not show_diff(out.getvalue(), """\
%#source info: s
%#  total number of:
%#    models:     2
%#    chains:     2
%#    alt. conf.: 1
%#    residues:   3
%#    atoms:      4
%#  number of atom element types: 1
%#  histogram of atom element frequency:
%#    "  " 4
%#  residue name classes:
%#    "common_amino_acid" 3
%#  number of chain ids: 1
%#  histogram of chain id frequency:
%#    "A" 2
%#  number of conformer ids: 3
%#  histogram of conformer id frequency:
%#    " " 1
%#    "A" 1
%#    "B" 1
%#  number of residue names: 1
%#  histogram of residue name frequency:
%#    "MET " 3
%#  number of groups of duplicate atom lables:    1
%#    total number of affected atoms:             2
%#    group " N   MET A   1 "
%#          " N   MET A   1 "
%#  hierarchy level of detail = atom
%#    model id=1 #chains=1
%#      chain id="A" #conformers=1
%#        conformer id=" " #residues=1 #atoms=2
%#          residue name="MET " seq=   1 icode=" " #atoms=2
%#             " N  "
%#             " N  "
%#    model id=2 #chains=1
%#      chain id="A" #conformers=2
%#        conformer id="A" #residues=1 #atoms=1 #alt. atoms=1
%#          residue name="MET " seq=   1 icode=" " #atoms=1
%#           & " N  "
%#        conformer id="B" #residues=1 #atoms=1 #alt. atoms=1
%#          residue name="MET " seq=   1 icode=" " #atoms=1
%#           & " N  "
""")
  #
  lines = flex.split_lines("""\
ATOM         CA  ASN     1
ATOM         C   ASN     1
ATOM         O   HOH     2
ATOM         H1  HOH     2
ATOM         H2  HOH     2
ATOM         CA  GLU     3
ATOM         C   GLU     3
ATOM         O   H2O     4
ATOM         O   OH2     5
ATOM         O   DOD     6
ATOM         P   U       7
ATOM         O   D2O     8
ATOM         O   OD2     9
ATOM         CA  ALA    10
ATOM         C   ALA    10
ATOM         O   TIP    11
ATOM         O   TIP3   12
""")
  hierarchy = pdb.input(source_info=None, lines=lines).construct_hierarchy()
  f = hierarchy.models()[0].chains()[0].conformers()[0]
  assert list(f.residue_class_selection(class_name="common_amino_acid")) == [0,2,9]
  assert list(f.residue_class_selection(class_name="common_rna_dna")) == [6]
  assert list(f.residue_class_selection(class_name="common_water")) \
      == [1,3,4,5,7,8,10,11]
  assert f.residue_class_selection(
    class_name="common_water", negate=False).all_eq(
      f.residue_class_selection(class_name="common_water"))
  assert list(f.residue_class_selection(
    class_name="common_water", negate=True)) == [0,2,6,9]
  assert list(f.residue_class_selection(residue_names=flex.std_string([]))) \
      == []
  assert list(f.residue_class_selection(residue_names=flex.std_string(["V"])))\
      == []
  assert list(f.residue_class_selection(
           residue_names=flex.std_string(["U", "OD2"]))) == [6,8]
  assert list(f.residue_class_selection(
                residue_names=flex.std_string(["U", "OD2"]), negate=True)) \
      == [0,1,2,3,4,5,7,9,10,11]
  #
  fw = f.select_residues(
    selection=f.residue_class_selection(class_name="common_water"))
  assert fw.residues_size() == 8
  atom_names = []
  for r in fw.residues():
    assert r.parent().memory_id() == fw.memory_id()
    assert r.atoms_size() > 0
    for a in r.atoms():
      assert a.parents_size() == 1
      assert a.parents()[0].memory_id() == r.memory_id()
      atom_names.append(a.name)
  assert ",".join(atom_names) \
      == " O  , H1 , H2 , O  , O  , O  , O  , O  , O  , O  "
  #
  for i in [0,1]:
    f.select_residues_in_place(
      selection=f.residue_class_selection(class_name="common_water"))
    out = StringIO()
    hierarchy.show(out=out)
    assert not show_diff(out.getvalue(), """\
model id=0 #chains=1
  chain id=" " #conformers=1
    conformer id=" " #residues=8 #atoms=10
      residue name="HOH " seq=   2 icode=" " #atoms=3
         " O  "
         " H1 "
         " H2 "
      residue name="H2O " seq=   4 icode=" " #atoms=1
         " O  "
      residue name="OH2 " seq=   5 icode=" " #atoms=1
         " O  "
      residue name="DOD " seq=   6 icode=" " #atoms=1
         " O  "
      residue name="D2O " seq=   8 icode=" " #atoms=1
         " O  "
      residue name="OD2 " seq=   9 icode=" " #atoms=1
         " O  "
      residue name="TIP " seq=  11 icode=" " #atoms=1
         " O  "
      residue name="TIP3" seq=  12 icode=" " #atoms=1
         " O  "
""")
  f.select_residues_in_place(selection=f.residue_class_selection(
    class_name="common_water", negate=True))
  out = StringIO()
  hierarchy.show(out=out)
  assert not show_diff(out.getvalue(), """\
model id=0 #chains=1
  chain id=" " #conformers=1
    conformer id=" " #residues=0 #atoms=0
""")
  hierarchy = pdb.input(source_info=None, lines=lines).construct_hierarchy()
  f = hierarchy.models()[0].chains()[0].conformers()[0]
  for i in [0,1]:
    f.select_residue_class_in_place(class_name="common_water", negate=True)
    out = StringIO()
    hierarchy.show(out=out)
    assert not show_diff(out.getvalue(), """\
model id=0 #chains=1
  chain id=" " #conformers=1
    conformer id=" " #residues=4 #atoms=7
      residue name="ASN " seq=   1 icode=" " #atoms=2
         " CA "
         " C  "
      residue name="GLU " seq=   3 icode=" " #atoms=2
         " CA "
         " C  "
      residue name="U   " seq=   7 icode=" " #atoms=1
         " P  "
      residue name="ALA " seq=  10 icode=" " #atoms=2
         " CA "
         " C  "
""")
  for i in [0,1]:
    f.select_residue_class_in_place(
      residue_names=flex.std_string(["GLU", "ALA"]), negate=False)
    out = StringIO()
    hierarchy.show(out=out)
    assert not show_diff(out.getvalue(), """\
model id=0 #chains=1
  chain id=" " #conformers=1
    conformer id=" " #residues=2 #atoms=4
      residue name="GLU " seq=   3 icode=" " #atoms=2
         " CA "
         " C  "
      residue name="ALA " seq=  10 icode=" " #atoms=2
         " CA "
         " C  "
""")
  for i in [0,1]:
    f.select_residue_class_in_place(
      residue_names=flex.std_string(["ALA"]), negate=True)
    out = StringIO()
    hierarchy.show(out=out)
    assert not show_diff(out.getvalue(), """\
model id=0 #chains=1
  chain id=" " #conformers=1
    conformer id=" " #residues=1 #atoms=2
      residue name="GLU " seq=   3 icode=" " #atoms=2
         " CA "
         " C  "
""")
  try: f.select_residue_class_in_place(class_name="xyz")
  except RuntimeError, e:
    assert not show_diff(str(e), 'unknown class_name="xyz"')
  else: raise RuntimeError("Exception expected.")
  try:
    f.select_residue_class_in_place(residue_names=flex.std_string(["vwxyz"]))
  except RuntimeError, e:
    assert not show_diff(str(e),
      'residue name with more than 4 characters: "vwxyz"')
  else: raise RuntimeError("Exception expected.")
  for s in [f.select_residues_in_place, f.select_residues]:
    try: s(flex.size_t([100]))
    except RuntimeError, e:
      assert not show_diff(str(e), "selection index out of range.")
    else: raise RuntimeError("Exception expected.")
  #
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM         CA  ASN     1
ATOM         P   +U      2
ATOM         O   HOH     3
ATOM        CD   CD      4
ATOM         S   SO4     5
ATOM         N   ABC     6
"""))
  check_hierarchy(
    hierarchy=pdb_inp.construct_hierarchy(),
    expected_formatted="""\
model id=0 #chains=1
  chain id=" " #conformers=1
    conformer id=" " #residues=6 #atoms=6
      residue name="ASN " seq=   1 icode=" " #atoms=1
         " CA "
      residue name="+U  " seq=   2 icode=" " #atoms=1
         " P  "
      residue name="HOH " seq=   3 icode=" " #atoms=1
         " O  "
      residue name="CD  " seq=   4 icode=" " #atoms=1
         "CD  "
      residue name="SO4 " seq=   5 icode=" " #atoms=1
         " S  "
      residue name="ABC " seq=   6 icode=" " #atoms=1
         " N  "
""",
    expected_overall_counts=dicts.easy(
      chain_ids={" ": 1},
      conformer_ids={" ": 1},
      residue_names={"ASN ": 1, "+U  ": 1, "HOH ": 1, "CD  ": 1,
                     "SO4 ": 1, "ABC ": 1},
      residue_name_classes={
        "common_water": 1,
        "common_element": 1,
        "common_rna_dna": 1,
        "common_amino_acid": 1,
        "common_small_molecule": 1,
        "other": 1}))
  #
  get_class = pdb.common_residue_names_get_class
  assert get_class(name="ALA") == "common_amino_acid"
  assert get_class(name="  U") == "common_rna_dna"
  assert get_class(name="HOH") == "common_water"
  assert get_class(name="SO4") == "common_small_molecule"
  assert get_class(name="CL ") == "common_element"
  assert get_class(name="ABC") == "other"

def exercise_xray_structure_simple():
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
CRYST1   61.410   54.829   43.543  90.00  90.00  90.00 P 21 21 21    8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.016284  0.000000  0.000000        0.00000
SCALE2      0.000000  0.018239  0.000000        0.00000
SCALE3      0.000000  0.000000  0.022966        0.00000
ATOM      1  N   GLN A   3      35.299  11.075  19.070  1.00 36.89           N
ATOM      2  CA  GLN A   3      34.482   9.927  18.794  0.63 37.88           C
ANISOU    2  CA  GLN A   3     7794   3221   3376  -1227   1064   2601       C
ATOM      3  Q   GLN A   3      35.130   8.880  17.864  0.84 37.52           C
ANISOU    3  Q   GLN A   3     7875   3041   3340   -981    727   2663       C
ATOM      4  O   GLN A   3      34.548   7.819  17.724  1.00 38.54      STUV
ATOM      5 1CB  GLN A   3      32.979  10.223  18.469  1.00 37.80
HETATM    6 CA   ION B   1      32.360  11.092  17.308  0.92 35.96          CA2+
HETATM    7 CA   ION B   2      30.822  10.665  17.190  1.00 36.87
"""))
  assert pdb_inp.atom_element_counts() == {" C": 2, "  ": 3, " N": 1, "CA": 1}
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
Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso
" N   GLN A   3 " N      4 ( 0.5748  0.2020  0.4380) 1.00 0.4672
" CA  GLN A   3 " C      4 ( 0.5615  0.1811  0.4316) 0.63 0.4797
     u_cart =  0.779  0.322  0.338 -0.123  0.106  0.260
" Q   GLN A   3 " C      4 ( 0.5721  0.1620  0.4103) 0.84 0.4752
     u_cart =  0.788  0.304  0.334 -0.098  0.073  0.266
" O   GLN A   3 " segid="STUV" O      4 ( 0.5626  0.1426  0.4070) 1.00 0.4881
"1CB  GLN A   3 " C      4 ( 0.5370  0.1865  0.4242) 1.00 0.4787
"CA   ION B   1 " Ca2+   4 ( 0.5270  0.2023  0.3975) 0.92 0.4554
"CA   ION B   2 " Ca     4 ( 0.5019  0.1945  0.3948) 1.00 0.4670
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
Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso
" N   GLN A   3 " N      1 (35.2990 11.0750 19.0700) 1.00 0.4672
" CA  GLN A   3 " C      1 (34.4820  9.9270 18.7940) 0.63 0.4797
     u_cart =  0.779  0.322  0.338 -0.123  0.106  0.260
" Q   GLN A   3 " C      1 (35.1300  8.8800 17.8640) 0.84 0.4752
     u_cart =  0.788  0.304  0.334 -0.098  0.073  0.266
" O   GLN A   3 " segid="STUV" O      1 (34.5480  7.8190 17.7240) 1.00 0.4881
"1CB  GLN A   3 " C      1 (32.9790 10.2230 18.4690) 1.00 0.4787
"CA   ION B   1 " Ca2+   1 (32.3600 11.0920 17.3080) 0.92 0.4554
"CA   ION B   2 " Ca     1 (30.8220 10.6650 17.1900) 1.00 0.4670
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
  except RuntimeError, e:
    assert not show_diff(str(e), """\
Unknown chemical element type: PDB ATOM " N   GLN A   3 " segid="STUV" element=" A"
  To resolve this problem, specify a chemical element type in
  columns 77-78 of the PDB file, right justified (e.g. " C").""")
  else: raise RuntimeError("Exception expected.")
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
  except RuntimeError, e:
    assert not show_diff(str(e), """\
Unknown chemical element type: PDB ATOM "1A   GLN A   3 " element="  "
  To resolve this problem, specify a chemical element type in
  columns 77-78 of the PDB file, right justified (e.g. " C").""")
  else: raise RuntimeError("Exception expected.")
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM      1  N   GLN A   3      35.299  11.075  99.070  1.00 36.89           B 5
"""))
  xray_structure = pdb_inp.xray_structure_simple(
    enable_scattering_type_unknown=True)
  assert xray_structure.scatterers()[0].scattering_type == "unknown"
  try:
    pdb_inp.xray_structure_simple()
  except RuntimeError, e:
    assert not show_diff(str(e), '''\
Unknown charge: PDB ATOM " N   GLN A   3 " element=" B" charge=" 5"''')
  else: raise RuntimeError("Exception expected.")
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
  except RuntimeError, e:
    assert not show_diff(str(e), '''\
Unknown scattering type: PDB ATOM " N   GLN A   3 " element="Cs" charge="3-"''')
  else: raise RuntimeError("Exception expected.")
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
  xray_structure = pdb_inp.xray_structure_simple()
  out = StringIO()
  xray_structure.show_scatterers(f=out)
  assert not show_diff(out.getvalue(), """\
Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso
" N   GLN A   3 " N      1 (35.2990 11.0750 19.0700) 1.00 0.4672
" CA  GLN A   3 " C      1 (34.4820  9.9270 18.7940) 0.63 0.4798
" N   GLN A   3 " N      1 (25.2990  1.0750  9.0700) 0.54 0.3406
" CA  GLN A   3 " C      1 (24.4820 -1.9270  8.7940) 1.00 0.3531
""")
  xray_structures = pdb_inp.xray_structures_simple()
  assert len(xray_structures) == 2
  out = StringIO()
  xray_structures[0].show_scatterers(f=out)
  assert not show_diff(out.getvalue(), """\
Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso
" N   GLN A   3 " N      1 (35.2990 11.0750 19.0700) 1.00 0.4672
" CA  GLN A   3 " C      1 (34.4820  9.9270 18.7940) 0.63 0.4798
""")
  out = StringIO()
  xray_structures[1].show_scatterers(f=out)
  assert not show_diff(out.getvalue(), """\
Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso
" N   GLN A   3 " N      1 (25.2990  1.0750  9.0700) 0.54 0.3406
" CA  GLN A   3 " C      1 (24.4820 -1.9270  8.7940) 1.00 0.3531
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
  pdb_inp = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
ATOM      1  N   GLN A   3      35.299  11.075  19.070  1.00 36.89           N
ATOM      2  CA  GLN A   3      34.482   9.927  18.794  0.63 37.88           C0
ATOM      3  C   GLN A   3      35.130   8.880  17.864  0.84 37.52           C 0
ATOM      4  O   GLN A   3      34.548   7.819  17.724  1.00 38.54           O00
ATOM      5 1CB  GLN A   3      32.979  10.223  18.469  1.00 37.80           C 1
HETATM    6 CA   ION B   1      32.360  11.092  17.308  0.92 35.96           X
HETATM    7 CA   ION B   2      30.822  10.665  17.190  1.00 36.87          FE4+
ATOM      8  O   MET A   5       6.215  22.789  24.067  1.00  0.00            -2
"""))
  assert [scatterer.scattering_type
    for scatterer in pdb_inp.xray_structure_simple(
      scattering_type_exact=True,
      enable_scattering_type_unknown=True).scatterers()] \
        == ["N", "C", "C", "O", "unknown", "unknown", "unknown", "O2-"]

def exercise(args):
  forever = "--forever" in args
  while True:
    exercise_atom()
    exercise_residue()
    exercise_chain()
    exercise_conformer()
    exercise_model()
    exercise_hierarchy()
    exercise_columns_73_76_evaluator()
    exercise_line_info_exceptions()
    exercise_pdb_input()
    exercise_xray_structure_simple()
    if (not forever): break
  print format_cpu_times()

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
