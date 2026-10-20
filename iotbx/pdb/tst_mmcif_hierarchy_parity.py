from __future__ import absolute_import, division, print_function
"""
Parity between the per-atom Python mmCIF hierarchy builder
(iotbx.pdb.mmcif.pdb_hierarchy_builder_python) and hierarchy_from_atom_site(),
the path the reader takes: compare_hierarchies() walks two live hierarchies in
lockstep (every stored attribute, every parent link, the flat atom order, the
derived labels), and the exercises below build both ways over the inline fixtures
of tst_mmcif_hierarchy.py, the model mmCIF files shipped with cctbx and the
phenix_regression models when present.
"""
import os
import copy
import glob
import libtbx.load_env
import iotbx.cif
from iotbx.pdb import mmcif
from cctbx.array_family import flex

LARGE_TIER_ATOMS = 500000

ATOM_ATTRIBUTES = ("serial", "name", "segid", "element", "charge", "hetero",
  "xyz", "sigxyz", "occ", "sigocc", "b", "sigb", "uij", "fp", "fdp",
  "resolution", "i_seq")

# extract_* arrays bound on hierarchy.atoms(); siguij and charge have no extractor
ARRAY_VEC = ("xyz", "sigxyz", "uij")
ARRAY_SCALAR = ("occ", "sigocc", "b", "sigb", "fp", "fdp", "resolution",
  "serial", "name", "segid", "element", "hetero", "i_seq")


class HierarchyMismatch(AssertionError):
  pass


def _check_eq(path, what, x, y):
  if x != y:
    raise HierarchyMismatch("%s: %s differs: %r != %r" % (path, what, x, y))


def _check_parent(path, child, parent):
  p = child.parent()
  if p is None or p.memory_id() != parent.memory_id():
    raise HierarchyMismatch("%s: parent link does not point at the enclosing node" % path)


def _compare_arrays(a, b):
  A = a.atoms()
  B = b.atoms()
  for what in ARRAY_VEC:
    x = getattr(A, "extract_" + what)().as_double()
    y = getattr(B, "extract_" + what)().as_double()
    if not x.all_eq(y):
      raise HierarchyMismatch("atoms().extract_%s() differs" % what)
  for what in ARRAY_SCALAR:
    x = getattr(A, "extract_" + what)()
    y = getattr(B, "extract_" + what)()
    if not x.all_eq(y):
      raise HierarchyMismatch("atoms().extract_%s() differs" % what)


def _compare_large(a, b):
  for im, (m1, m2) in enumerate(zip(a.models(), b.models())):
    pm = "model[%d] id=%r" % (im, m1.id)
    _check_eq(pm, "id", m1.id, m2.id)
    if not m1.is_identical_hierarchy(m2):
      raise HierarchyMismatch("%s: is_identical_hierarchy is False" % pm)
  _compare_arrays(a, b)
  # root.residue_groups() is a generator (Python injector); model and chain
  # accessors return lists
  r1 = list(a.residue_groups())
  r2 = list(b.residue_groups())
  _check_eq("root", "residue_groups_size", len(r1), len(r2))
  for i, (x, y) in enumerate(zip(r1, r2)):
    for what in ("resseq", "icode", "link_to_previous"):
      _check_eq("residue_groups()[%d]" % i, what, getattr(x, what), getattr(y, what))


def compare_hierarchies(a, b, tier=None):
  """
  Raise HierarchyMismatch at the first difference between the two hierarchies.
  tier "full": lockstep walk over models, chains, residue groups, atom groups and
  atoms comparing every stored attribute and every parent back-link, then the flat
  atoms() order against the walk, atoms_with_labels() element by element, and the
  extract_* arrays. tier "large": per-model is_identical_hierarchy, the extract_*
  arrays, and a loop over residue_groups(). Default tier by atom count.
  """
  n = a.atoms_size()
  _check_eq("root", "atoms_size", n, b.atoms_size())
  if tier is None:
    tier = "large" if n > LARGE_TIER_ATOMS else "full"
  ma = a.models()
  mb = b.models()
  _check_eq("root", "models_size", len(ma), len(mb))
  if tier == "large":
    _compare_large(a, b)
    return
  walked = []
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
          pg = "%s ag[%d] altloc=%r resname=%r" % (pr, ig, ag1.altloc, ag1.resname)
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
            walked.append((x.memory_id(), y.memory_id()))
  for i, (x, y) in enumerate(zip(a.atoms(), b.atoms())):
    if x.memory_id() != walked[i][0] or y.memory_id() != walked[i][1]:
      raise HierarchyMismatch("atoms()[%d] is not the %d-th atom of the tree walk" % (i, i))
  for i, (x, y) in enumerate(zip(a.atoms_with_labels(), b.atoms_with_labels())):
    p = "atoms_with_labels()[%d]" % i
    _check_eq(p, "id_str", x.id_str(), y.id_str())
    _check_eq(p, "is_first_in_chain", x.is_first_in_chain, y.is_first_in_chain)
    _check_eq(p, "is_first_after_break", x.is_first_after_break, y.is_first_after_break)
  _compare_arrays(a, b)


def _first_block(text):
  model = iotbx.cif.reader(input_string=text).model()
  return list(model.values())[0]


two_residue_block = """\
data_t
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
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
  block = _first_block(two_residue_block)
  a = mmcif.pdb_hierarchy_builder_python(block).hierarchy
  b = mmcif.pdb_hierarchy_builder_python(block).hierarchy
  compare_hierarchies(a, b)                 # identical builds pass
  compare_hierarchies(a, b, tier="large")   # the large tier passes too
  # a numeric difference is found with the path to the atom
  b.atoms()[3].b = 13.5
  try:
    compare_hierarchies(a, b)
  except HierarchyMismatch as e:
    assert "atom[0] name=' N  '" in str(e) and "b differs" in str(e), str(e)
  else:
    raise AssertionError("mismatch not detected")
  try:
    compare_hierarchies(a, b, tier="large")
  except HierarchyMismatch as e:
    assert "extract_b" in str(e), str(e)
  else:
    raise AssertionError("mismatch not detected on the large tier")
  b.atoms()[3].b = 13.0
  compare_hierarchies(a, b)
  # a wrong parent link is found even though every attribute agrees: the two CA
  # atoms (identical fields, altloc A and B) are swapped between their groups
  # without re-parenting, so counts, names and values match the original
  c = mmcif.pdb_hierarchy_builder_python(block).hierarchy
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
  except HierarchyMismatch as e:
    assert "parent link" in str(e), str(e)
  else:
    raise AssertionError("broken parent link not detected")
  # link_to_previous is compared on both tiers
  d = mmcif.pdb_hierarchy_builder_python(block).hierarchy
  d.models()[0].chains()[0].residue_groups()[1].link_to_previous = False
  for tier in ("full", "large"):
    try:
      compare_hierarchies(a, d, tier=tier)
    except HierarchyMismatch as e:
      assert "link_to_previous" in str(e), str(e)
    else:
      raise AssertionError("link_to_previous difference not detected on tier %s" % tier)


def _sorted_copy(h):
  """The three calls cif_input.construct_hierarchy makes after the build."""
  s = h.deep_copy()
  s.sort_atoms_in_place()
  s.reset_atom_i_seqs()
  s.atoms_reset_serial()
  return s


def _with_ins_code_column(block):
  """
  Parity twin for inputs without _atom_site.pdbx_PDB_ins_code: the Python
  builder gives icode '' there (None reaches the residue_group constructor),
  but ' ' for files that carry the column with '?' or '.'. The parity loop
  compares on a copy with the column added as '?', so that a builder giving
  ' ' in both cases compares equal. Returns the block itself when the column
  is present.
  """
  if "_atom_site.pdbx_PDB_ins_code" in block:
    return block
  block = copy.deepcopy(block)   # the caller's block stays unchanged
  loop = block.get_loop_or_row("_atom_site")
  loop.add_column("_atom_site.pdbx_PDB_ins_code",
    flex.std_string(loop.size(), "?"))
  return block


def compare_block(block, label, tier=None):
  """
  Build block with both builders and compare, unsorted and sorted. If the
  reference raises, the new path must raise too. Returns "ok", "both-raise" or
  "f6" (compared on the insertion-code twin).
  """
  outcome = "ok"
  if "_atom_site.pdbx_PDB_ins_code" not in block:
    block = _with_ins_code_column(block)
    outcome = "f6"
  try:
    ref = mmcif.pdb_hierarchy_builder_python(block).hierarchy
  except Exception as e_ref:
    try:
      mmcif.hierarchy_from_atom_site(block)
    except Exception:
      return "both-raise"
    raise AssertionError("%s: reference raised %s: %s, new path did not" % (
      label, type(e_ref).__name__, e_ref))
  new = mmcif.hierarchy_from_atom_site(block)
  try:
    compare_hierarchies(ref, new, tier=tier)
    compare_hierarchies(_sorted_copy(ref), _sorted_copy(new), tier=tier)
  except HierarchyMismatch as e:
    raise HierarchyMismatch("%s: %s" % (label, e))
  return outcome


# Expected compare_block outcome of every inline fixture. The three "f6" ones
# have no insertion-code column; input_missing_mandatory_items lacks mandatory
# items, so both builders must raise. Anything else is "ok". Pinned so that a
# fixture that starts raising on both sides is not mistaken for a comparison.
EXPECTED_FIXTURE_OUTCOME = {
  "input_missing_mandatory_items": "both-raise",
  "input_1ab1": "f6",
  "input_4edr": "f6",
  "input_null_occ_b": "f6",
}


def parity_fixtures():
  """
  Every inline mmCIF model fixture of tst_mmcif_hierarchy.py, as (name, text).
  A fixture that is only an _atom_site loop body is prefixed with
  pdb_atom_site_loop_header, as its exercise there does.
  """
  from iotbx.pdb import tst_mmcif_hierarchy as t
  result = []
  for name in ("input_1ab1", "input_4edr", "input_1ezu", "input_charges",
               "input_anom", "input_multi_model_single_chain",
               "input_question_mark_altloc", "input_null_occ_b",
               "input_missing_mandatory_items"):
    result.append((name, getattr(t, name)))
  for name in ("input_4ehz", "input_3zdi", "input_4gln", "input_1ezu_loop",
               "input_2hok"):
    result.append((name, t.pdb_atom_site_loop_header + getattr(t, name)))
  return result


def exercise_fixtures():
  n = 0
  for name, text in parity_fixtures():
    outcome = compare_block(_first_block(text), "fixture " + name)
    expected = EXPECTED_FIXTURE_OUTCOME.get(name, "ok")
    assert outcome == expected, "fixture %s: outcome %r, expected %r" % (
      name, outcome, expected)
    n += 1
  assert n >= 12, n


def _repo_model_files():
  result = []
  for module, rel in (("iotbx", "regression/data"), ("mmtbx", "regression/pdbs"),
                      ("xcif", "regression")):
    d = libtbx.env.dist_path(module)
    if d is None:
      continue
    for fn in sorted(glob.glob(os.path.join(d, rel, "*.cif"))):
      if os.path.basename(fn) == "example.cif":
        # no pdbx_PDB_model_num: pdb_hierarchy_builder_python raises TypeError
        continue
      with open(fn, "rb") as f:   # bytes: a stray non-UTF-8 byte must not abort the scan
        if b"_atom_site.Cartn_x" in f.read():
          result.append(fn)
  return result


def compare_file(fn, tier=None):
  model = iotbx.cif.reader(file_path=fn).model()
  return compare_block(list(model.values())[0], os.path.basename(fn), tier=tier)


def exercise_repo_files():
  files = _repo_model_files()
  assert len(files) >= 9, files
  for fn in files:
    outcome = compare_file(fn)
    assert outcome != "both-raise", "%s: both builders raised, nothing compared" % fn


def exercise_phenix_regression():
  d = libtbx.env.find_in_repositories("phenix_regression")
  if d is None:
    print("phenix_regression not present, skipping")
    return
  n = 0
  outcomes = {}
  for root, dirs, files in os.walk(d):
    for name in sorted(files):
      if not name.endswith(".cif"):
        continue
      fn = os.path.join(root, name)
      with open(fn, "rb") as f:
        if b"_atom_site.Cartn_x" not in f.read():
          continue
      outcome = compare_file(fn)
      outcomes[outcome] = outcomes.get(outcome, 0) + 1
      n += 1
  print("phenix_regression: %d model mmCIF files compared, outcomes %s" % (n, outcomes))
  assert n >= 300, n
  assert outcomes.get("both-raise", 0) == 0, (
    "files on which both builders raised were not compared: %s" % outcomes)


def run():
  exercise_compare_hierarchies()
  exercise_fixtures()
  exercise_repo_files()
  exercise_phenix_regression()


if __name__ == "__main__":
  run()
  print("OK")
