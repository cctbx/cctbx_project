"""
Parity between the per-atom Python mmCIF hierarchy builder
(iotbx.pdb.mmcif._pdb_hierarchy_builder_python) and hierarchy_from_atom_site(),
the path the reader takes: compare_hierarchies() walks two live hierarchies in
lockstep (every stored attribute, every parent link, the extract_* arrays),
and the exercises below build both ways over the inline fixtures of
tst_mmcif_hierarchy.py and a fixed list of model mmCIF files from the
repositories (PARITY_FILES). The reader's sort is not compared: it reads only
residue and atom names, which the walk compares.
"""
import os
import libtbx.load_env
from libtbx.test_utils import Exception_expected
import iotbx.cif
from iotbx.pdb import mmcif
from iotbx.pdb import tst_mmcif_hierarchy as fixtures

ATOM_ATTRIBUTES = ("serial", "name", "segid", "element", "charge", "hetero",
  "xyz", "sigxyz", "occ", "sigocc", "b", "sigb", "uij", "fp", "fdp",
  "resolution", "i_seq")

# extract_* arrays bound on hierarchy.atoms(); siguij and charge have no
# extractor
ARRAY_VEC = ("xyz", "sigxyz", "uij")
ARRAY_SCALAR = ("occ", "sigocc", "b", "sigb", "fp", "fdp", "resolution",
  "serial", "name", "segid", "element", "hetero", "i_seq")


def _eq(x, y):
  """x == y, with NaN equal to NaN (xyz, sigxyz and uij are tuples)."""
  return x == y or repr(x) == repr(y)


def _check_eq(path, what, x, y):
  if not _eq(x, y):
    raise AssertionError("%s: %s differs: %r != %r" % (path, what, x, y))


def _check_parent(path, child, parent):
  p = child.parent()
  if p is None or p.memory_id() != parent.memory_id():
    raise AssertionError(
      "%s: parent link does not point at the enclosing node" % path)


def _arrays_eq(x, y):
  """x.all_eq(y), with NaN equal to NaN."""
  return x.all_eq(y) or repr(list(x)) == repr(list(y))


def _compare_arrays(a, b):
  A = a.atoms()
  B = b.atoms()
  for what in ARRAY_VEC:
    x = getattr(A, "extract_" + what)().as_double()
    y = getattr(B, "extract_" + what)().as_double()
    if not _arrays_eq(x, y):
      raise AssertionError("atoms().extract_%s() differs" % what)
  for what in ARRAY_SCALAR:
    x = getattr(A, "extract_" + what)()
    y = getattr(B, "extract_" + what)()
    if not _arrays_eq(x, y):
      raise AssertionError("atoms().extract_%s() differs" % what)


def compare_hierarchies(a, b):
  """
  Raise AssertionError at the first difference between two hierarchies:
  a lockstep walk over models, chains, residue groups, atom groups and atoms
  comparing every stored attribute and every parent back-link, then the
  extract_* arrays.

  Parameters
  ----------
  a, b : iotbx.pdb.hierarchy.root
    The hierarchies, a the reference.

  Raises
  ------
  AssertionError
    With the path to the first node that differs.
  """
  _check_eq("root", "atoms_size", a.atoms_size(), b.atoms_size())
  ma = a.models()
  mb = b.models()
  _check_eq("root", "models_size", len(ma), len(mb))
  for im, (m1, m2) in enumerate(zip(ma, mb)):
    pm = "model[%d] id=%r" % (im, m1.id)
    _check_eq(pm, "id", m1.id, m2.id)
    _check_parent(pm, m1, a)
    _check_parent(pm, m2, b)
    c1 = m1.chains()
    c2 = m2.chains()
    _check_eq(pm, "chains_size", len(c1), len(c2))
    for ic, (ch1, ch2) in enumerate(zip(c1, c2)):
      pc = "%s chain[%d] id=%r" % (pm, ic, ch1.id)
      _check_eq(pc, "id", ch1.id, ch2.id)
      _check_parent(pc, ch1, m1)
      _check_parent(pc, ch2, m2)
      r1 = ch1.residue_groups()
      r2 = ch2.residue_groups()
      _check_eq(pc, "residue_groups_size", len(r1), len(r2))
      for ir, (rg1, rg2) in enumerate(zip(r1, r2)):
        pr = "%s rg[%d] resid=%r" % (pc, ir, rg1.resid())
        for what in ("resseq", "icode", "link_to_previous"):
          _check_eq(pr, what, getattr(rg1, what), getattr(rg2, what))
        _check_parent(pr, rg1, ch1)
        _check_parent(pr, rg2, ch2)
        g1 = rg1.atom_groups()
        g2 = rg2.atom_groups()
        _check_eq(pr, "atom_groups_size", len(g1), len(g2))
        for ig, (ag1, ag2) in enumerate(zip(g1, g2)):
          pg = "%s ag[%d] altloc=%r resname=%r" % (
            pr, ig, ag1.altloc, ag1.resname)
          _check_eq(pg, "altloc", ag1.altloc, ag2.altloc)
          _check_eq(pg, "resname", ag1.resname, ag2.resname)
          _check_parent(pg, ag1, rg1)
          _check_parent(pg, ag2, rg2)
          at1 = ag1.atoms()
          at2 = ag2.atoms()
          _check_eq(pg, "atoms_size", len(at1), len(at2))
          for ia, (x, y) in enumerate(zip(at1, at2)):
            pa = "%s atom[%d] name=%r" % (pg, ia, x.name)
            for what in ATOM_ATTRIBUTES:
              _check_eq(pa, what, getattr(x, what), getattr(y, what))
            _check_parent(pa, x, ag1)
            _check_parent(pa, y, ag2)
  _compare_arrays(a, b)


two_residue_block = fixtures.atom_site_header + """\
ATOM 1 N N . ALA A 1 5 ? 1.0 2.0 3.0 1.00 10.00 5 ALA A N 1
ATOM 2 C CA A ALA A 1 5 ? 1.5 2.0 3.0 0.50 11.00 5 ALA A CA 1
ATOM 2 C CA B ALA A 1 5 ? 1.5 2.0 3.0 0.50 11.00 5 ALA A CA 1
ATOM 4 N N . GLY A 1 6 ? 2.0 2.0 3.0 1.00 13.00 6 GLY A N 1
HETATM 5 O O . HOH B 2 . ? 9.0 9.0 9.0 1.00 30.00 101 HOH A O 1
"""
# rows 2 and 3 differ only in altloc on purpose (same id, so the serial agrees
# too): the parent-link case below swaps them between their atom groups, so
# every attribute still agrees and only the links differ.


def exercise_compare_hierarchies():
  block = fixtures.first_block(two_residue_block)
  a = mmcif._pdb_hierarchy_builder_python(block).hierarchy
  b = mmcif._pdb_hierarchy_builder_python(block).hierarchy
  compare_hierarchies(a, b)                 # identical builds pass
  # a numeric difference is found with the path to the atom
  b.atoms()[3].b = 13.5
  try:
    compare_hierarchies(a, b)
  except AssertionError as e:
    assert "atom[0] name=' N  '" in str(e) and "b differs" in str(e), str(e)
  else:
    raise Exception_expected
  b.atoms()[3].b = 13.0
  compare_hierarchies(a, b)
  # NaN equals NaN, in an attribute and in an array; NaN against a number is
  # a difference
  nan = float("nan")
  a.atoms()[3].b = nan
  b.atoms()[3].b = nan
  compare_hierarchies(a, b)
  b.atoms()[3].b = 13.0
  try:
    compare_hierarchies(a, b)
  except AssertionError as e:
    assert "b" in str(e), str(e)
  else:
    raise Exception_expected
  a.atoms()[3].b = 13.0
  compare_hierarchies(a, b)
  # a wrong parent link is found even though every attribute agrees: the two CA
  # atoms (identical fields, altloc A and B) are swapped between their groups
  # without re-parenting, so counts, names and values match the original
  c = mmcif._pdb_hierarchy_builder_python(block).hierarchy
  rg = c.models()[0].chains()[0].residue_groups()[0]
  ag_a, ag_b = rg.atom_groups()[1], rg.atom_groups()[2]
  atom_a, atom_b = ag_a.atoms()[0], ag_b.atoms()[0]
  ag_a.remove_atom(i=0)
  ag_b.remove_atom(i=0)
  ag_a.append_atom_with_other_parent(atom_b)
  ag_b.append_atom_with_other_parent(atom_a)
  assert c.atoms_size() == a.atoms_size()
  try:
    compare_hierarchies(a, c)
  except AssertionError as e:
    assert "parent link" in str(e), str(e)
  else:
    raise Exception_expected
  # link_to_previous is compared
  d = mmcif._pdb_hierarchy_builder_python(block).hierarchy
  d.models()[0].chains()[0].residue_groups()[1].link_to_previous = False
  try:
    compare_hierarchies(a, d)
  except AssertionError as e:
    assert "link_to_previous" in str(e), str(e)
  else:
    raise Exception_expected


def compare_block(block, label):
  """Build block with both builders and compare; a mismatch names label."""
  ref = mmcif._pdb_hierarchy_builder_python(block).hierarchy
  new = mmcif.hierarchy_from_atom_site(block)
  try:
    compare_hierarchies(ref, new)
  except AssertionError as e:
    raise AssertionError("%s: %s" % (label, e))


def parity_fixtures():
  """
  The inline mmCIF model fixtures of tst_mmcif_hierarchy.py that both builders
  build, as (name, text). A fixture that is only an _atom_site loop body is
  prefixed with its header, as its exercise there does.
  """
  result = []
  for name in ("input_1ab1", "input_4edr", "input_1ezu", "input_charges",
               "input_anom", "input_multi_model_single_chain",
               "input_question_mark_altloc", "input_null_occ_b"):
    result.append((name, getattr(fixtures, name)))
  for name in ("input_4ehz", "input_3zdi", "input_4gln", "input_1ezu_loop",
               "input_2hok"):
    result.append((name,
      fixtures.pdb_atom_site_loop_header + getattr(fixtures, name)))
  result.append(("input_two_atoms",
    fixtures.atom_site_header + fixtures.input_two_atoms))
  return result


def exercise_fixtures():
  for name, text in parity_fixtures():
    compare_block(fixtures.first_block(text), "fixture " + name)


# The model files compared, as (module, path within the module), each for
# what it exercises. The parity_*.cif files hold the residue groups that
# carry a feature (first altloc, insertion code, hybrid-36 residue number,
# two-letter element, HETATM, anisotrop row) and the first residue groups of
# every chain, with their _atom_site_anisotrop rows, cell and symmetry. With
# the fixtures the test takes about a second (about 4 s with PYTHONDEVMODE=1
# and PYTHONTRACEMALLOC=1, as the unix CI jobs run it). Formal charges, B
# tensors, fp/fdp and null occupancies are covered by the fixtures.
PARITY_FILES = [
  # PDB entry: label_asym_id chain split after the polymer, HETATM
  ("xcif", "regression/1yjp.cif"),
  # no model number, no insertion-code column, label columns only
  ("xcif", "regression/example.cif"),
  # predicted model: label columns only, one chain
  ("iotbx", "regression/data/non_zero_origin_ncs_model.cif"),
  # deposited entry: anisotrop U with esd columns, two-letter elements,
  # chain split before the ligand and the waters
  ("iotbx", "regression/data/parity_3orl.cif"),
  # insertion codes, altlocs, hybrid-36 residue numbers, segids
  ("iotbx", "regression/data/parity_icode_altloc.cif"),
  # two models with the same atom ids, altlocs
  ("iotbx", "regression/data/parity_two_models.cif"),
  # protein chain split from its ligand, hybrid-36 residue numbers
  ("iotbx", "regression/data/parity_protein_ligand.cif"),
  # a ligand alone, HETATM
  ("iotbx", "regression/data/parity_ligand_only.cif"),
  # a CA trace: one atom per residue group, three chains
  ("iotbx", "regression/data/parity_ca_trace.cif"),
]


def compare_file(fn):
  """compare_block on the first data block of the file, as the reader does."""
  model = iotbx.cif.reader(file_path=fn).model()
  compare_block(list(model.values())[0], os.path.basename(fn))


def exercise_model_files():
  for module, rel in PARITY_FILES:
    fn = os.path.join(libtbx.env.dist_path(module), rel)
    assert os.path.isfile(fn), fn
    compare_file(fn)


def run():
  exercise_compare_hierarchies()
  exercise_fixtures()
  exercise_model_files()


if __name__ == "__main__":
  run()
  print("OK")
