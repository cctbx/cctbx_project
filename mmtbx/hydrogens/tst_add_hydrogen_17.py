"""place_hydrogens with exclude_water=False and a water without O: its H are
removed and reported as missing a neighbouring heavy atom; waters with O are
unchanged. workaround_002() places water H only for a residue with exactly
one O and two H/D, never from another residue's O."""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from mmtbx.hydrogens import reduce_hydrogen
from six.moves import cStringIO as StringIO

cryst1 = "CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1\n"

def water(resseq, x, altloc=" ", o=True, o_altloc=" "):
  lines = []
  if o:
    lines.append("HETATM    1  O  %sHOH A%4d    %8.3f  10.000  10.000  1.00 20.00"
                 "           O" % (o_altloc, resseq, x))
  for name, dx, dy in (("H1", 0.9, 0.1), ("H2", -0.3, 0.9)):
    lines.append("HETATM    1  %-3s%sHOH A%4d    %8.3f%8.3f  10.000  1.00 20.00"
                 "           H" % (name, altloc, resseq, x + dx, 10.0 + dy))
  return "\n".join(lines) + "\n"

def place(pdb_str, keep_existing_H=True):
  model = mmtbx.model.manager(model_input=iotbx.pdb.input(
    lines=(cryst1 + pdb_str).split("\n"), source_info=None), log=null_out())
  obj = reduce_hydrogen.place_hydrogens(model=model, exclude_water=False,
    keep_existing_H=keep_existing_H)
  obj.run()
  log = StringIO()
  obj.show(log, verbose=True)
  return obj, obj.get_model().get_hierarchy(), log.getvalue()

def water_atoms(hierarchy, resseq):
  return [a for a in hierarchy.atoms()
          if a.parent().parent().resseq.strip() == str(resseq)]

def check_ideal(hierarchy, resseq, d=0.85):
  """per conformer: O and two H at d, 103.91 deg"""
  atoms = water_atoms(hierarchy, resseq)
  o = [a for a in atoms if a.element.strip() == "O"]
  assert len(o) == 1, resseq
  hs = [a for a in atoms if a.element.strip() == "H"]
  assert len(hs) % 2 == 0 and len(hs) > 0, (resseq, len(hs))
  O = matrix.col(o[0].xyz)
  by_alt = {}
  for a in hs:
    by_alt.setdefault(a.parent().altloc, []).append(matrix.col(a.xyz) - O)
  for alt, vs in by_alt.items():
    assert len(vs) == 2, (resseq, alt)
    for v in vs:
      assert approx_equal(v.length(), d, eps=1e-3, out=None), v.length()
    assert approx_equal(vs[0].angle(vs[1], deg=True), 103.91, eps=0.01,
                        out=None)
  return sorted((a.name, a.parent().altloc, a.xyz) for a in hs)

def check_reported(obj, report, resseq):
  labels = [" H1  HOH A%4d " % resseq, " H2  HOH A%4d " % resseq]
  for label in labels:
    assert label in obj.site_labels_missing_neighbor, (
      label, obj.site_labels_missing_neighbor)
    assert label not in obj.site_labels_no_para
  assert obj.residues_missing_neighbor.count("HOH A %d" % resseq) == 1
  assert "HOH A %d" % resseq in report.split(
    "were not placed because a neighbouring heavy atom is")[1]
  assert "could not be parameterized" not in report

def test_first_water_without_o():
  # (a) the H of a water without O come before any O: no crash, removed
  obj, h, report = place(water(1, 5.0, o=False) + water(2, 15.0))
  assert water_atoms(h, 1) == []
  check_ideal(h, 2)
  check_reported(obj, report, 1)

def test_water_without_o_after_normal():
  # (b) after a normal water: removed, and the normal water as on its own
  obj, h, report = place(water(1, 5.0) + water(2, 15.0, o=False))
  assert water_atoms(h, 2) == []
  check_reported(obj, report, 2)
  alone = check_ideal(place(water(1, 5.0))[1], 1)
  assert check_ideal(h, 1) == alone

def test_normal_waters():
  # (c) waters with O are unchanged by the fix. A normal water gets ideal H
  # around its own O. A blank O with H1/H2 in altloc A (existing behaviour,
  # kept): placement adds a blank H1/H2 pair, workaround_002() then sees four
  # H in conformer A and leaves the A pair at its input positions.
  obj, h, report = place(water(1, 5.0) + water(2, 15.0, altloc="A"))
  check_ideal(h, 1)
  w2 = dict(((a.name.strip(), a.parent().altloc), a.xyz)
            for a in water_atoms(h, 2) if a.element.strip() == "H")
  assert sorted(w2) == [("H1", ""), ("H1", "A"), ("H2", ""), ("H2", "A")]
  assert approx_equal(w2[("H1", "A")], (15.9, 10.1, 10.0))
  assert approx_equal(w2[("H2", "A")], (14.7, 10.9, 10.0))
  o = matrix.col([a for a in water_atoms(h, 2) if a.name.strip() == "O"][0].xyz)
  for name in ("H1", "H2"):
    d = (matrix.col(w2[(name, "")]) - o).length()
    assert approx_equal(d, 0.85, eps=1e-3, out=None), d
  # the A pair has no restraint next to the blank pair of the same names; as
  # before, it is listed as not parameterized but kept
  assert obj.site_labels_missing_neighbor == []
  assert obj.site_labels_no_para == [" H1 AHOH A   2 ", " H2 AHOH A   2 "]

def test_keep_existing_h_false():
  # (d) keep_existing_H=False: input H are replaced; a water without O has
  # nothing left, waters with O get ideal H
  for pdb_str, with_o, without_o in (
      (water(1, 5.0, o=False) + water(2, 15.0), [2], [1]),
      (water(1, 5.0) + water(2, 15.0, o=False), [1], [2]),
      (water(1, 5.0) + water(2, 15.0, altloc="A"), [1, 2], [])):
    obj, h, report = place(pdb_str, keep_existing_H=False)
    for resseq in with_o:
      check_ideal(h, resseq)
    for resseq in without_o:
      assert water_atoms(h, resseq) == []

if __name__ == "__main__":
  test_first_water_without_o()
  test_water_without_o_after_normal()
  test_normal_waters()
  test_keep_existing_h_false()
  print("OK")
