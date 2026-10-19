"""place_hydrogens with exclude_water=False and a water without O: its H are
removed and reported as missing a neighbouring heavy atom; waters with O are
unchanged. workaround_002() places water H only for a residue with exactly
one O and two H/D, never from another residue's O.

A residue with a blank atom_group and exactly one altloc is one conformer:
missing H are added once, to the blank group, and existing altloc H are kept
(no duplicate names). Blank + two altlocs: H go into each altloc."""
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
  # around its own O. A blank O with H1/H2 in altloc A: the A pair is kept
  # and idealized, no blank pair is added
  obj, h, report = place(water(1, 5.0) + water(2, 15.0, altloc="A"))
  check_ideal(h, 1)
  w2 = check_ideal(h, 2)
  assert [(name, alt) for name, alt, xyz in w2] == [
    (" H1 ", "A"), (" H2 ", "A")], w2
  assert obj.site_labels_missing_neighbor == []
  assert obj.site_labels_no_para == []

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

# Ala-Ser-Lys; Ser OG is moved into altlocs by the tests
pdb_peptide = """\
ATOM      1  N   ALA A   1      11.204   2.315  10.312  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.555   2.835  10.528  1.00 20.00           C
ATOM      3  C   ALA A   1      13.588   1.726  10.381  1.00 20.00           C
ATOM      4  O   ALA A   1      13.325   0.567  10.721  1.00 20.00           O
ATOM      5  CB  ALA A   1      12.866   3.965   9.556  1.00 20.00           C
ATOM      6  N   SER A   2      14.767   2.090   9.884  1.00 20.00           N
ATOM      7  CA  SER A   2      15.853   1.130   9.701  1.00 20.00           C
ATOM      8  C   SER A   2      16.968   1.370  10.713  1.00 20.00           C
ATOM      9  O   SER A   2      17.281   2.516  11.041  1.00 20.00           O
ATOM     10  CB  SER A   2      16.411   1.224   8.278  1.00 20.00           C
ATOM     11  OG  SER A   2      15.456   0.826   7.311  1.00 20.00           O
ATOM     12  N   LYS A   3      17.560   0.282  11.199  1.00 20.00           N
ATOM     13  CA  LYS A   3      18.643   0.370  12.172  1.00 20.00           C
ATOM     14  C   LYS A   3      19.945  -0.148  11.570  1.00 20.00           C
ATOM     15  O   LYS A   3      19.964  -1.176  10.892  1.00 20.00           O
ATOM     16  CB  LYS A   3      18.290  -0.424  13.432  1.00 20.00           C
ATOM     17  OXT LYS A   3      20.995   0.451  11.763  1.00 20.00           O
"""
og_line = "ATOM     11  OG  SER A   2      15.456   0.826   7.311  1.00"

def ser_og_in_a(occ=0.6):
  """blank Ser backbone and CB, OG only in altloc A"""
  assert pdb_peptide.count(og_line) == 1
  return pdb_peptide.replace(og_line,
    "ATOM     11  OG ASER A   2      15.456   0.826   7.311  %.2f" % occ)

def ser_og_in_a_b():
  """blank Ser backbone and CB, OG in altlocs A and B"""
  return ser_og_in_a(occ=0.5).replace("ATOM     12  N   LYS",
    "ATOM     11  OG BSER A   2      16.300   0.000   7.300  0.50 20.00"
    "           O\nATOM     12  N   LYS")

def with_h(pdb_str):
  """the model with H placed from scratch, as PDB records"""
  return place(pdb_str, keep_existing_H=False)[1].as_pdb_string()

def ser_h(hierarchy):
  """{altloc: sorted H names} of Ser 2"""
  result = {}
  for ag in hierarchy.atom_groups():
    if ag.resname != "SER": continue
    result[ag.altloc] = sorted(a.name.strip() for a in ag.atoms()
                               if a.element.strip() == "H")
  return result

def ser_atom(hierarchy, name, altloc):
  atoms = [a for a in hierarchy.atoms() if a.parent().resname == "SER" and
           a.name.strip() == name and a.parent().altloc == altloc]
  assert len(atoms) == 1, (name, altloc, len(atoms))
  return atoms[0]

def check_hg(hierarchy, hg_altloc, og_altloc):
  hg = ser_atom(hierarchy, "HG", hg_altloc)
  og = ser_atom(hierarchy, "OG", og_altloc)
  d = (matrix.col(hg.xyz) - matrix.col(og.xyz)).length()
  assert approx_equal(d, 0.84, eps=1e-3, out=None), d
  return hg

def test_single_altloc_side_chain_kept():
  # (c) Ser with OG and HG only in altloc A (0.6), keep_existing_H=True: the
  # input HG stays in A with its OG, no second HG; backbone and CB H in the
  # blank group
  # H from the model without altlocs, then OG and HG moved to altloc A
  lines = []
  hg_xyz = None
  for line in with_h(pdb_peptide).splitlines():
    if line[12:16] in (" OG ", " HG ") and line[17:20] == "SER":
      line = line[:16] + "A" + line[17:54] + "  0.60" + line[60:]
      if line[12:16] == " HG ":
        hg_xyz = tuple(float(line[i:i+8]) for i in (30, 38, 46))
    lines.append(line)
  pdb_str = "\n".join(lines) + "\n"
  assert pdb_str.count("ASER") == 2 and hg_xyz is not None
  obj, h, report = place(pdb_str)
  assert ser_h(h) == {"": ["H", "HA", "HB2", "HB3"], "A": ["HG"]}, ser_h(h)
  hg = check_hg(h, "A", "A")
  assert approx_equal(hg.xyz, hg_xyz, eps=1e-3)
  assert approx_equal(hg.occ, 0.6)
  assert obj.site_labels_no_para == []

def test_single_altloc_side_chain_placed():
  # (e) the same Ser without input H, keep_existing_H=False: as before the
  # fix. All new H go into the blank group, also HG whose OG is in A.
  # Placing each new H in its parent's atom_group would put HG into A
  # (open point "Where placed H go in a single-altloc residue")
  obj, h, report = place(ser_og_in_a(), keep_existing_H=False)
  assert ser_h(h) == {"": ["H", "HA", "HB2", "HB3", "HG"], "A": []}, ser_h(h)
  hg = check_hg(h, "", "A")
  assert approx_equal(hg.occ, 0.6)
  assert obj.site_labels_no_para == []

def test_two_altlocs_unchanged():
  # (b) water with a blank O and H1/H2 in A and in B, (d) Ser with OG in A
  # and B: blank-parent H go into each altloc, as before the fix
  pdb_str = """\
HETATM    1  O   HOH A   1       5.000  10.000  10.000  1.00 20.00           O
HETATM    2  H1 AHOH A   1       5.900  10.100  10.000  0.50 20.00           H
HETATM    3  H2 AHOH A   1       4.700  10.900  10.000  0.50 20.00           H
HETATM    4  H1 BHOH A   1       4.100  10.000  10.000  0.50 20.00           H
HETATM    5  H2 BHOH A   1       5.000   9.100  10.000  0.50 20.00           H
"""
  obj, h, report = place(pdb_str)
  assert [(name, alt) for name, alt, xyz in check_ideal(h, 1)] == [
    (" H1 ", "A"), (" H1 ", "B"), (" H2 ", "A"), (" H2 ", "B")]
  assert obj.site_labels_no_para == []
  obj, h, report = place(pdb_str, keep_existing_H=False)
  assert [(name, alt) for name, alt, xyz in check_ideal(h, 1)] == [
    (" H1 ", ""), (" H2 ", "")]
  heavy = ser_og_in_a_b()
  h_input = with_h(heavy)
  no_hg = "\n".join(l for l in h_input.splitlines()
                    if not (l[12:16] == " HG " and l[17:20] == "SER"))
  assert h_input.count(" HG ") == 2 and no_hg.count(" HG ") == 0
  conformer_h = ["H", "HA", "HB2", "HB3", "HG"]
  for pdb_str, keep in ((h_input, True), (no_hg, True), (heavy, False)):
    obj, h, report = place(pdb_str, keep_existing_H=keep)
    assert ser_h(h) == {"": [], "A": conformer_h, "B": conformer_h}, ser_h(h)
    check_hg(h, "A", "A")
    check_hg(h, "B", "B")
    assert obj.site_labels_no_para == []

if __name__ == "__main__":
  test_first_water_without_o()
  test_water_without_o_after_normal()
  test_normal_waters()
  test_keep_existing_h_false()
  test_single_altloc_side_chain_kept()
  test_single_altloc_side_chain_placed()
  test_two_altlocs_unchanged()
  print("OK")
