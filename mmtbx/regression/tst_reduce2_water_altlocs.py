"""reduce2 qmready=True labels unlabelled partial-occupancy water alternates
(iotbx.pdb.water_altlocs, model.label_water_alternates) before placing H:
two waters 0.8 A apart at occupancy 0.5 are conformers A and B, each with
two H, instead of the "Bad clash(es)" Sorry; labels only (the non-H atoms
keep their occupancies and B; the qmready regularization may move them). Without qmready nothing is labelled."""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
from scitbx import matrix
from mmtbx.regression.tst_reduce2_capping import run, check_sorry, pdb_lys_asp_gly

# ASP A 5 - GLY A 6 of pdb_lys_asp_gly (complete) and two waters near OD2 of
# ASP 5 (12.531 8.382 6.067), 0.8 A apart, at occupancy 0.5 each
waters = """\
HETATM   20  O   HOH A 101      14.900   9.900   5.000  0.50 20.00           O
HETATM   21  O   HOH A 102      15.700   9.900   5.000  0.50 20.00           O
"""
pdb_str = "".join(l + "\n" for l in pdb_lys_asp_gly.splitlines()
  if not (l.startswith("ATOM") and l[17:20] == "LYS")).replace("END\n",
  waters + "END\n")

def heavy(hierarchy):
  return dict(((a.parent().parent().resseq.strip(), a.name.strip()),
    (a.parent().altloc, a.xyz, a.occ, a.b)) for a in hierarchy.atoms()
    if not a.element_is_hydrogen())

def water_h(hierarchy, resseq):
  return [(a.parent().altloc, a.name.strip()) for a in hierarchy.atoms()
    if a.parent().parent().resseq.strip() == resseq and a.element_is_hydrogen()]

def test_qmready_labels():
  prefix = "tst_reduce2_water_altlocs_1"
  fn_out = prefix + "_out.pdb"
  r = run(pdb_str, prefix, ["qmready=True"], fn_out=fn_out)
  h = iotbx.pdb.input(fn_out).construct_hierarchy()
  before = heavy(iotbx.pdb.input(source_info=None,
    lines=pdb_str.splitlines()).construct_hierarchy())
  after = heavy(h)
  assert after[("101", "O")][0] == "A" and after[("102", "O")][0] == "B", \
    (after[("101", "O")], after[("102", "O")])
  # labels only: occupancies and B of every non-H atom as read, sites
  # within what the qmready regularization moves
  for key, (altloc, xyz, occ, b) in before.items():
    if key[0] in ("101", "102"): continue
    assert after[key][0] == altloc, key
  for key in before:
    candidates = [before[key]]
    if key[1] in ("OD1", "OD2"):
      # reduce2 may swap the carboxylate O names (IUPAC convention)
      candidates = [before[(key[0], n)] for n in ("OD1", "OD2")]
    altloc, xyz, occ, b = min(candidates, key=lambda c:
      (matrix.col(after[key][1]) - matrix.col(c[1])).length())
    assert (matrix.col(after[key][1]) - matrix.col(xyz)).length() < 0.3, key
    assert after[key][2] == occ and after[key][3] == b, key
  # each conformer's water has two H
  assert sorted(water_h(h, "101")) == [("A", "H1"), ("A", "H2")], water_h(h, "101")
  assert sorted(water_h(h, "102")) == [("B", "H1"), ("B", "H2")], water_h(h, "102")
  assert r.total_charge == -1, r.total_charge
  with open(prefix + "_out.txt") as f:
    description = f.read()
  assert "Water alternates" in description and "HOH A 101 -> altloc A" in \
    description and "HOH A 102 -> altloc B" in description, description

def test_without_qmready_unchanged():
  prefix = "tst_reduce2_water_altlocs_2"
  fn_out = prefix + "_out.pdb"
  run(pdb_str, prefix, [], fn_out=fn_out)
  after = heavy(iotbx.pdb.input(fn_out).construct_hierarchy())
  assert after[("101", "O")][0].strip() == "" and \
    after[("102", "O")][0].strip() == "", after
  with open(prefix + "_out.txt") as f:
    assert "Water alternates" not in f.read()

def test_full_occupancy_still_a_clash():
  # two waters at full occupancy on overlapping sites are an error, not
  # alternates: the clash check stops as before
  text = pdb_str.replace("  0.50 20.00           O", "  1.00 20.00           O")
  check_sorry(text, "tst_reduce2_water_altlocs_3", ["qmready=True"],
    "Bad clash(es) prevents qmready=True", "HOH A 101", "HOH A 102")

if __name__ == "__main__":
  test_qmready_labels()
  test_without_qmready_unchanged()
  test_full_occupancy_still_a_clash()
  print("OK")
