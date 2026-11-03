"""Test iotbx.pdb.altloc_completion and root.complete_altlocs"""
import iotbx.pdb
from iotbx.pdb import altloc_completion
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry
from six.moves import cStringIO as StringIO

def atom_line(name, altloc, resname, resseq, xyz, element, chain="A",
              occ=0.50):
  return "ATOM      1 %-4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f 20.00" \
    "          %2s" % (" "+name if len(name) < 4 else name, altloc, resname,
    chain, resseq, xyz[0], xyz[1], xyz[2], occ, element)

def hierarchy(lines):
  return iotbx.pdb.input(source_info=None, lines=lines).construct_hierarchy()

def xxx(atoms, resname="XXX"):
  """atoms: list of (name, altloc, xyz, occ), all carbon"""
  return [atom_line(name, altloc, resname, 1, xyz, "C", occ=occ)
    for name, altloc, xyz, occ in atoms]

def atom_named(h, name, altloc=None):
  result = [a for a in h.atoms() if a.name.strip() == name and
    (altloc is None or a.parent().altloc == altloc)]
  assert len(result) == 1, result
  return result[0]

def occupancies_of(h, name):
  return dict([(a.parent().altloc, a.occ) for a in h.atoms()
    if a.name.strip() == name])

def single_source_lines(resname="XXX"):
  return xxx([
    ("C1", "",  (0,   0,   0), 1.00),
    ("C2", "A", (1.5, 0,   0), 0.46),
    ("C3", "A", (3.0, 0,   0), 0.46),
    ("C3", "B", (3.0, 1.5, 0), 0.54)], resname=resname)

def multi_source_lines(occ_c2_a=0.81, occ_c3_d=0.20):
  # 1R6J SER 205 pattern
  return xxx([
    ("C1", "",  (0,   0,   0), 1.00),
    ("C2", "A", (1.5, 0,   0), occ_c2_a),
    ("C2", "B", (1.5, 0.7, 0), 0.19),
    ("C3", "A", (2.9, 0,   0), 0.43),
    ("C3", "B", (2.9, 0.7, 0), 0.10),
    ("C3", "C", (2.2, 1.1, 0), 0.18),
    ("C3", "D", (2.5, 1.3, 0), occ_c3_d)])

def exercise_single_source_shared():
  h = hierarchy(single_source_lines())
  rg_id = " XXX A   1 "
  r = h.complete_altlocs()
  assert r.shared == [(rg_id, "C2", "A")]
  assert r.copied == []
  assert r.truncated == []
  assert r.changed
  a = atom_named(h, "C2")
  assert a.parent().altloc == ""
  assert approx_equal(a.occ, 0.46)
  assert altloc_completion.report_lines(r) == [
    "%s C2: altloc A -> shared (needed by B)" % rg_id]
  assert list(h.atoms().extract_i_seq()) == list(range(h.atoms_size()))

def exercise_multi_source_copies():
  h = hierarchy(multi_source_lines())
  rg_id = " XXX A   1 "
  r = h.complete_altlocs()
  copied = sorted(r.copied)
  assert [c[:4] for c in copied] == [
    (rg_id, "C2", "A", "C"), (rg_id, "C2", "A", "D")]
  assert approx_equal([c[4] for c in copied], [0.18, 0.20])
  occ = occupancies_of(h, "C2")
  assert sorted(occ.keys()) == ["A", "B", "C", "D"]
  assert approx_equal([occ[k] for k in "ABCD"], [0.43, 0.19, 0.18, 0.20])
  assert approx_equal(atom_named(h, "C2", "C").xyz, (1.5, 0, 0))
  assert r.shared == []
  assert r.inconsistent == []
  assert r.changed
  assert altloc_completion.report_lines(r) == [
    "%s C2: copied A -> C" % rg_id, "%s C2: copied A -> D" % rg_id]

def exercise_inconsistent_occupancies():
  # occ_R of A, C, D sum to 0.81, source has 0.90: split proportionally
  h = hierarchy(multi_source_lines(occ_c2_a=0.90))
  rg_id = " XXX A   1 "
  r = h.complete_altlocs()
  assert r.inconsistent == [rg_id]
  occ = occupancies_of(h, "C2")
  assert approx_equal(occ["A"] + occ["C"] + occ["D"], 0.90)
  assert approx_equal([occ[k] for k in "ACD"],
    [0.90*0.43/0.81, 0.90*0.18/0.81, 0.90*0.20/0.81])

def exercise_ambiguous_sorry():
  h = hierarchy(multi_source_lines(occ_c3_d=0.45))
  s0 = h.as_pdb_string()
  try: h.complete_altlocs()
  except Sorry as e:
    s = str(e)
    assert 'residue " XXX A   1 "' in s, s
    assert "atom C2 could come from" in s and "conformer D" in s, s
    assert "altloc A (1.50 A from C1, occ 0.81)" in s, s
    assert "altloc B (1.17 A from C3, occ 0.19)" in s, s
    assert "fix the model by hand" in s, s
  else: raise AssertionError("Sorry expected")
  assert h.as_pdb_string() == s0

def exercise_shared_then_copy_into_same_group():
  # CD is label C's only atom and is shared (C atom group empties), then
  # CG A is copied into C: the copy must not be lost
  lines = xxx([
    ("C1", "",  (0,   0, 0), 1.0),
    ("CG", "A", (3,   0, 0), 0.4),
    ("CG", "B", (0,   3, 0), 0.4),
    ("CD", "C", (1.5, 0, 0), 0.2)])
  h = hierarchy(lines)
  totals0 = {}
  for a in h.atoms():
    totals0[a.name.strip()] = totals0.get(a.name.strip(), 0) + a.occ
  r = h.complete_altlocs()
  assert r.shared == [(" XXX A   1 ", "CD", "C")]
  assert [c[:4] for c in r.copied] == [(" XXX A   1 ", "CG", "A", "C")]
  assert h.atoms_size() == 5
  totals1 = {}
  for a in h.atoms():
    totals1[a.name.strip()] = totals1.get(a.name.strip(), 0) + a.occ
  assert sorted(totals1.keys()) == sorted(totals0.keys())
  for name in totals0:
    assert approx_equal(totals1[name], totals0[name])
  occ = occupancies_of(h, "CG")
  assert approx_equal([occ[k] for k in "ABC"], [0.4*2/3, 0.4, 0.4/3])
  assert atom_named(h, "CD").parent().altloc == ""
  assert [ag.altloc for ag in h.only_residue_group().atom_groups()] == \
    ["", "A", "B", "C"]

def exercise_truncated():
  # 1R6J ASP 195 pattern
  h = hierarchy(xxx([
    ("C1", "A", (0,   0, 0), 0.5),
    ("C1", "B", (0,   3, 0), 0.5),
    ("C2", "A", (1.5, 0, 0), 0.5)]))
  rg_id = " XXX A   1 "
  s0 = h.as_pdb_string()
  r = h.complete_altlocs()
  assert r.truncated == [(rg_id, "B", ["C2"], "C2: A copy is 3.4 A from C1 B")]
  assert r.shared == [] and r.copied == [] and not r.changed
  assert h.as_pdb_string() == s0
  assert altloc_completion.report_lines(r) == [
    "%s C2: A copy is 3.4 A from C1 B, not borrowed; truncated in B" % rg_id]

def exercise_capping_h_not_borrowed():
  # H11 B caps C1 B where B is truncated (no C2): bonded to C1 A, but its
  # parent C1 B is not in A, so A does not borrow it
  lines = xxx([
    ("C1", "A", (0,   0,   0), 0.5),
    ("C1", "B", (0,   2.2, 0), 0.5),
    ("C2", "A", (1.5, 0,   0), 0.5)])
  lines.append(atom_line("H11", "B", "XXX", 1, (0, 1.2, 0), "H"))
  h = hierarchy(lines)
  assert altloc_completion.bonded(atom_named(h, "H11"), atom_named(h, "C1",
    "A"))
  s0 = h.as_pdb_string()
  r = h.complete_altlocs()
  assert r.shared == [] and r.copied == [] and not r.changed
  # the gap H11 stays open in A and is reported like any open gap
  rg_id = " XXX A   1 "
  assert [t[:3] for t in r.truncated] == [
    (rg_id, "A", ["H11"]), (rg_id, "B", ["C2"])]
  assert h.as_pdb_string() == s0

def exercise_cap_h_is_no_anchor():
  # H11 B caps C1 B (1.0 A) and is 1.0 A from C2 A: C2 A is bonded in B
  # only through H11, so B does not borrow it, nor A H11
  lines = xxx([
    ("C1", "A", (0,   0,   0), 0.5),
    ("C1", "B", (1.5, 2.3, 0), 0.5),
    ("C2", "A", (1.5, 0,   0), 0.5)])
  lines.append(atom_line("H11", "B", "XXX", 1, (1.5, 1.3, 0), "H"))
  h = hierarchy(lines)
  h11, c2a = atom_named(h, "H11"), atom_named(h, "C2")
  assert altloc_completion.bonded(h11, c2a)
  assert not altloc_completion.bonded(atom_named(h, "C1", "B"), c2a)
  s0 = h.as_pdb_string()
  r = h.complete_altlocs()
  assert r.shared == [] and r.copied == [] and not r.changed
  rg_id = " XXX A   1 "
  assert [t[:3] for t in r.truncated] == [
    (rg_id, "A", ["H11"]), (rg_id, "B", ["C2"])]
  assert h.as_pdb_string() == s0

def exercise_zero_occupancy_sorry():
  # C borrows C2 from A and A borrows C3 from C (copies: B has C2 only)
  h = hierarchy(xxx([
    ("C2", "A", (0,   0, 0), 0),
    ("C2", "B", (0,   3, 0), 0),
    ("C3", "C", (1.5, 0, 0), 0)]))
  s0 = h.as_pdb_string()
  try: h.complete_altlocs()
  except Sorry as e:
    assert str(e) == 'all conformers of residue " XXX A   1 " have zero ' \
      'occupancy; cannot split', str(e)
  else: raise AssertionError("Sorry expected")
  assert h.as_pdb_string() == s0

def exercise_hydrogens_follow():
  lines = single_source_lines()
  lines.insert(2, atom_line("H21", "A", "XXX", 1, (1.5, 1.0, 0), "H",
    occ=0.46))
  lines.insert(3, atom_line("H22", "A", "XXX", 1, (1.5, -1.0, 0), "H",
    occ=0.46))
  h = hierarchy(lines)
  r = h.complete_altlocs()
  assert sorted(n for _, n, _ in r.shared) == ["C2", "H21", "H22"]
  for name in ["C2", "H21", "H22"]:
    assert atom_named(h, name).parent().altloc == ""
  assert r.copied == [] and r.truncated == []

def exercise_hydrogen_only_label():
  lines = xxx([
    ("C1", "",  (0,   0,   0), 1.0),
    ("C2", "A", (1.5, 0,   0), 0.5),
    ("C2", "B", (1.5, 1.0, 0), 0.5)])
  lines.append(atom_line("H11", "C", "XXX", 1, (-0.5, 0.9, 0), "H"))
  h = hierarchy(lines)
  rg = h.only_residue_group()
  assert altloc_completion.conformer_labels(rg) == ["A", "B"]
  s0 = h.as_pdb_string()
  r = h.complete_altlocs()
  assert not r.changed and r.truncated == []
  assert h.as_pdb_string() == s0

def exercise_single_label():
  h = hierarchy(xxx([
    ("C1", "",  (0,   0, 0), 1.0),
    ("C2", "A", (1.5, 0, 0), 0.5),
    ("C3", "A", (3.0, 0, 0), 0.5)]))
  r = h.complete_altlocs()
  assert not r.changed and r.truncated == []

def exercise_water_and_ligand():
  h = hierarchy([
    atom_line("O", "A", "HOH", 1, (0, 0, 0), "O"),
    atom_line("O", "B", "HOH", 1, (0, 1.2, 0), "O")])
  r = h.complete_altlocs()
  assert not r.changed and r.truncated == []
  h = hierarchy(single_source_lines(resname="LIG"))
  rg_id = " LIG A   1 "
  r = h.complete_altlocs()
  assert r.shared == [(rg_id, "C2", "A")]
  assert atom_named(h, "C2").parent().altloc == ""

def exercise_microheterogeneity():
  h = hierarchy([
    atom_line("N",  "",  "PRO", 1, (0,   0,   0), "N", occ=1.0),
    atom_line("CA", "A", "PRO", 1, (1.4, 0,   0), "C"),
    atom_line("C",  "A", "PRO", 1, (2.0, 1.3, 0), "C"),
    atom_line("CB", "A", "PRO", 1, (2.0, -1.3, 0), "C"),
    atom_line("CA", "B", "SER", 1, (1.4, 0.3, 0), "C"),
    atom_line("C",  "B", "SER", 1, (2.0, 1.6, 0), "C")])
  s0 = h.as_pdb_string()
  r = h.complete_altlocs()
  assert not r.changed and r.truncated == []
  assert h.as_pdb_string() == s0

def exercise_check_only_and_idempotent():
  h = hierarchy(single_source_lines())
  s0 = h.as_pdb_string()
  r1 = h.complete_altlocs(check_only=True)
  assert r1.changed
  assert h.as_pdb_string() == s0
  log = StringIO()
  r2 = h.complete_altlocs(log=log)
  assert r2.shared == r1.shared
  assert log.getvalue().splitlines() == altloc_completion.report_lines(r2)
  assert h.as_pdb_string() != s0
  r3 = h.complete_altlocs()
  assert not r3.changed and r3.shared == [] and r3.copied == []
  assert altloc_completion.report_lines(r3) == []

def exercise_scattering_neutral():
  lines = ["CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1"]
  lines += multi_source_lines()
  def f_calc(xrs):
    return xrs.structure_factors(d_min=1.5, algorithm="direct").f_calc()
  pdb_inp = iotbx.pdb.input(source_info=None, lines=lines)
  f0 = f_calc(pdb_inp.xray_structure_simple())
  h = pdb_inp.construct_hierarchy()
  r = h.complete_altlocs()
  assert len(r.copied) == 2
  xrs = h.extract_xray_structure(
    crystal_symmetry=pdb_inp.crystal_symmetry())
  assert xrs.scatterers().size() == 9
  f1 = f_calc(xrs)
  assert f1.indices().all_eq(f0.indices())
  from cctbx.array_family import flex
  d = flex.abs(f1.data() - f0.data())
  assert flex.max(d) < 1.e-6 * flex.max(flex.abs(f0.data()))

def exercise_label_occupancy_and_bonded():
  h = hierarchy(multi_source_lines())
  rg = h.only_residue_group()
  assert approx_equal(altloc_completion.label_occupancy(rg, "A"), 0.43)
  assert approx_equal(altloc_completion.label_occupancy(rg, "B"), 0.10)
  assert altloc_completion.conformer_labels(rg) == ["A", "B", "C", "D"]
  h = hierarchy(xxx([
    ("C1", "", (0,    0, 0), 1.0),
    ("C2", "", (1.9,  0, 0), 1.0),
    ("C3", "", (0, 1.93, 0), 1.0)]))
  c1, c2, c3 = [atom_named(h, n) for n in ["C1", "C2", "C3"]]
  assert altloc_completion.bonded(c1, c2)
  assert not altloc_completion.bonded(c1, c3)

def run():
  exercise_single_source_shared()
  exercise_multi_source_copies()
  exercise_inconsistent_occupancies()
  exercise_ambiguous_sorry()
  exercise_shared_then_copy_into_same_group()
  exercise_truncated()
  exercise_capping_h_not_borrowed()
  exercise_cap_h_is_no_anchor()
  exercise_zero_occupancy_sorry()
  exercise_hydrogens_follow()
  exercise_hydrogen_only_label()
  exercise_single_label()
  exercise_water_and_ligand()
  exercise_microheterogeneity()
  exercise_check_only_and_idempotent()
  exercise_scattering_neutral()
  exercise_label_occupancy_and_bonded()
  print("OK")

if (__name__ == "__main__"):
  run()
