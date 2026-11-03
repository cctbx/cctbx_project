"""Test residue_group.relabel_altlocs and root.reconcile_altlocs"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
from cctbx import crystal
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from six.moves import cStringIO as StringIO

def atom_line(name, altloc, resname, resseq, xyz, element, chain="A"):
  return "ATOM      1 %-4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f  0.50 20.00" \
    "          %2s" % (" "+name if len(name) < 4 else name, altloc, resname,
    chain, resseq, xyz[0], xyz[1], xyz[2], element)

def water(resseq, sites):
  """sites: list of (altloc, xyz)"""
  return [atom_line("O", altloc, "HOH", resseq, xyz, "O")
    for altloc, xyz in sites]

def hierarchy(lines):
  return iotbx.pdb.input(source_info=None, lines=lines).construct_hierarchy()

def altloc_at(h, resseq, xyz):
  """altloc of the atom of residue resseq located at xyz"""
  for atom in h.atoms():
    if (atom.parent().parent().resseq_as_int() == resseq and
        approx_equal(atom.xyz, xyz, eps=1.e-3, out=None)):
      return atom.parent().altloc
  raise RuntimeError("atom not found")

def exercise_relabel_altlocs():
  h = hierarchy(water(1, [("A", (0,0,0)), ("B", (6,0,0)), ("C", (0,6,0))]))
  rg = h.only_residue_group()
  rg.relabel_altlocs({"A":"B", "B":"A"})
  assert altloc_at(h, 1, (0,0,0)) == "B"
  assert altloc_at(h, 1, (6,0,0)) == "A"
  assert altloc_at(h, 1, (0,6,0)) == "C"
  # atom order is kept (arrays aligned with the hierarchy stay valid)
  assert approx_equal(h.atoms().extract_xyz(), [(0,0,0), (6,0,0), (0,6,0)])
  rg.relabel_altlocs({"A":"C", "B":"A", "C":"B"})
  assert approx_equal(h.atoms().extract_xyz(), [(0,0,0), (6,0,0), (0,6,0)])
  assert altloc_at(h, 1, (0,0,0)) == "A"
  assert altloc_at(h, 1, (6,0,0)) == "C"
  assert altloc_at(h, 1, (0,6,0)) == "B"
  # mapping must be a permutation of non-blank labels
  for bad in [{"A":"B"}, {"A":"C", "C":"C"}, {"":"A", "A":""}]:
    try: rg.relabel_altlocs(bad)
    except AssertionError: pass
    else: raise RuntimeError("Exception expected: %s" % str(bad))
  # rejected mapping leaves labels untouched
  assert altloc_at(h, 1, (0,0,0)) == "A"
  # labels absent from the residue_group are allowed
  rg.relabel_altlocs({"D":"E", "E":"D"})
  assert altloc_at(h, 1, (6,0,0)) == "C"

def two_inconsistent_waters():
  # Water 1 A clashes with water 2 A, water 1 B with water 2 B (2 A apart).
  # Swapping labels of either water removes both clashes.
  return water(1, [("A", (0,0,0)), ("B", (6,0,0))]) + \
         water(2, [("A", (0,2,0)), ("B", (6,2,0))])

def exercise_two_waters_swapped():
  h = hierarchy(two_inconsistent_waters())
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 1
  assert len(r.unresolved) == 0
  # swap the later residue, keep the first one
  assert r.swaps[0][0].strip() == "A   2"
  assert r.swaps[0][1] == {"A":"B", "B":"A"}
  assert altloc_at(h, 1, (0,0,0)) == "A"
  assert altloc_at(h, 2, (6,2,0)) == "A"
  assert altloc_at(h, 2, (0,2,0)) == "B"
  # idempotent: second run changes nothing
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 0
  assert len(r.unresolved) == 0

def exercise_check_only():
  h = hierarchy(two_inconsistent_waters())
  s = h.as_pdb_string()
  log = StringIO()
  r = h.reconcile_altlocs(crystal_symmetry=None, check_only=True, log=log)
  assert len(r.swaps) == 1
  assert h.as_pdb_string() == s
  assert log.getvalue().find("A   2") >= 0

def exercise_distant_sites_untouched():
  # Same as two_inconsistent_waters but water 2 moved 10 A away
  h = hierarchy(water(1, [("A", (0,0,0)), ("B", (6,0,0))]) +
                water(2, [("A", (0,12,0)), ("B", (6,12,0))]))
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert r.n_sites == 2
  assert len(r.swaps) == 0
  assert altloc_at(h, 2, (0,12,0)) == "A"

def exercise_fewest_changes():
  # Chain of contacts 1-2-3; water 2 is inconsistent with both neighbours.
  # Swapping water 2 alone is better than swapping waters 1 and 3.
  h = hierarchy(
    water(1, [("A", (0,0,0)), ("B", (6,0,0))]) +
    water(2, [("A", (0,2,0)), ("B", (6,2,0))]) +
    water(3, [("A", (0,4,0)), ("B", (6,4,0))]))
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert [s[0].strip() for s in r.swaps] == ["A   2"]
  assert altloc_at(h, 2, (0,2,0)) == "B"
  assert altloc_at(h, 2, (6,2,0)) == "A"

def exercise_three_labels():
  # Water 2 atom P clashes with water 1 A and B, atom Q with water 1 B,
  # atom R with nothing: the only consistent labels are P=C, Q=A, R=B.
  P, Q, R = (1.8,1,0), (5.6,0,0), (0,20,0)
  h = hierarchy(
    water(1, [("A", (0,0,0)), ("B", (3.6,0,0)), ("C", (0,10,0))]) +
    water(2, [("A", P), ("B", Q), ("C", R)]))
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 1
  assert len(r.unresolved) == 0
  assert altloc_at(h, 2, P) == "C"
  assert altloc_at(h, 2, Q) == "A"
  assert altloc_at(h, 2, R) == "B"

def exercise_unresolvable_flagged():
  # Three waters, every pair inconsistent (frustrated triangle): any
  # relabeling leaves a clash. Must be reported and left unchanged.
  h = hierarchy(
    water(1, [("A", (0,0,0)),   ("B", (20,0,0))]) +
    water(2, [("A", (2,0,0)),   ("B", (22,0,0))]) +
    water(3, [("A", (1,1.7,0)), ("B", (21,1.7,0))]))
  s = h.as_pdb_string()
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 0
  assert len(r.unresolved) == 1
  assert sorted([i.strip() for i in r.unresolved[0]]) == \
    ["A   1", "A   2", "A   3"]
  assert h.as_pdb_string() == s

def exercise_unresolvable_kept_local():
  # Frustrated triangle (waters 1-3, A and B) linked through water 4 to an
  # inconsistent pair (waters 5 and 6) far away: the pair is still fixed.
  h = hierarchy(
    water(1, [("A", (0,0,0)),   ("B", (20,0,0))]) +
    water(2, [("A", (2,0,0)),   ("B", (22,0,0))]) +
    water(3, [("A", (1,1.7,0)), ("B", (21,1.7,0))]) +
    water(4, [("A", (1,3.7,0)), ("B", (30,0,0))]) +
    water(5, [("A", (30,2,0)),  ("B", (50,0,0))]) +
    water(6, [("A", (50,2,0)),  ("B", (60,0,0))]))
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.unresolved) == 1
  assert ["A   1", "A   2", "A   3"] == \
    sorted([i.strip() for i in r.unresolved[0]])[:3]
  assert "A   6" in [i.strip() for i, m in r.swaps]
  assert altloc_at(h, 6, (50,2,0)) == "B"
  for resseq in (1, 2, 3):
    assert altloc_at(h, resseq, (20+2*(resseq==2)+(resseq==3),
      1.7*(resseq==3), 0)) == "B"

def exercise_marginal_frustration_ignored():
  # Frustrated triangle with marginal overlaps (2.3 A): not reported
  h = hierarchy(
    water(1, [("A", (0,0,0)),       ("B", (20,0,0))]) +
    water(2, [("A", (2.3,0,0)),     ("B", (22.3,0,0))]) +
    water(3, [("A", (1.15,1.992,0)), ("B", (21.15,1.992,0))]))
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.unresolved) == 0

def exercise_marginal_clash_ignored():
  # Water 2 is marginally too close (2.38 A) to water 1 A whatever its
  # labels are: a tight contact, not a labeling problem.
  h = hierarchy(water(1, [("A", (0,0,0)), ("B", (10,0,0))]) +
                water(2, [("A", (2.38,0,0)), ("B", (0,2.38,0))]))
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 0
  assert len(r.unresolved) == 0

def exercise_insignificant_gain_ignored():
  # Tight contact between nearly identical conformers: swapping water 2
  # lowers the overlap by only 0.002 A. Not worth a relabeling.
  h = hierarchy(water(1, [("A", (0,0,0)), ("B", (-0.1,0,0))]) +
                water(2, [("A", (2.3,0,0)), ("B", (2.3,0.4,0))]))
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 0
  assert len(r.unresolved) == 0

def exercise_symmetry_contact():
  # Water 2 A clashes with water 1 A only through a lattice translation
  cs = crystal.symmetry((10,10,10,90,90,90), "P 1")
  lines = water(1, [("A", (1,5,5)), ("B", (5,5,5))]) + \
          water(2, [("A", (9,5,5)), ("B", (5,8.5,5))])
  h = hierarchy(lines)
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 0
  h = hierarchy(lines)
  r = h.reconcile_altlocs(crystal_symmetry=cs, log=null_out())
  assert len(r.swaps) == 1
  assert altloc_at(h, 2, (9,5,5)) == "B"

def exercise_nonsense_symmetry():
  # Placeholder or incomplete symmetry (e.g. CRYST1 1 1 1 in cryo-EM models)
  # is treated as an isolated molecule
  for cs in [crystal.symmetry((1,1,1,90,90,90), "P 1"),
             crystal.symmetry(space_group_symbol="P 1"),
             crystal.symmetry(unit_cell=(50,50,50,90,90,90))]:
    h = hierarchy(two_inconsistent_waters())
    r = h.reconcile_altlocs(crystal_symmetry=cs, log=null_out())
    assert len(r.swaps) == 1
    assert len(r.unresolved) == 0
    assert altloc_at(h, 2, (6,2,0)) == "A"

def exercise_clash_through_other_conformer():
  # C of residue 1 A clashes with N of residue 2 A. Both are bonded to C of
  # residue 1 B (peptide link B-A), but B is another conformer, so the pair
  # is not 1-3 excluded: the clash is scored.
  from iotbx.pdb import altloc_reconciliation
  h = hierarchy(
    [atom_line("C", "A", "ALA", 1, (0,0,0), "C"),
     atom_line("C", "B", "ALA", 1, (1,1,0), "C"),
     atom_line("N", "A", "ALA", 2, (2,0,0), "N"),
     atom_line("O1", "B", "ALA", 2, (-1.9,0,0), "O"),
     atom_line("O2", "B", "ALA", 2, (0,-1.9,0), "O"),
     atom_line("O3", "B", "ALA", 2, (-1.3,-1.3,0), "O")])
  atoms = altloc_reconciliation._atoms(h.only_model(), False)
  edges = altloc_reconciliation._edge_terms(atoms=atoms,
    symmetry=altloc_reconciliation._symmetry(None, atoms.sites_cart, 4.5),
    distance_cutoff=4.5, clash_tolerance=0.5)
  assert edges[(0,1)][("A","A")][0] > 0.5
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 0

def exercise_geometry_clash_not_flagged():
  # Water 2 clashes with water 1 under any labels (bad geometry): relabeling
  # cannot help, so this is not a labeling problem.
  h = hierarchy(water(1, [("A", (0,0,0)), ("B", (0.3,0,0))]) +
                water(2, [("A", (2.3,0,0)), ("B", (2.6,0,0))]))
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 0
  assert len(r.unresolved) == 0

def exercise_partial_labels_flagged():
  # Waters 2 and 3 have only A; they clash with water 1 A and B
  # respectively. Water 2 would fit as B, but labels can only be permuted
  # within a residue: reported.
  h = hierarchy(water(1, [("A", (0,0,0)), ("B", (10,0,0))]) +
                water(2, [("A", (2,0,0))]) +
                water(3, [("A", (12,0,0))]))
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 0
  assert len(r.unresolved) == 1

def dipeptide(swap_second):
  # ALA-ALA, both residues with A and B conformers; B = A shifted by 4 A.
  # Wrong pairing breaks the C(1)-N(2) peptide bond without any clash.
  res1 = [("N", (-1.521, 3.898, 0.713)), ("CA", (-1.879, 2.678, 1.370)),
          ("CB", (-3.346, 2.715, 1.747)), ("C", (-1.551, 1.486, 0.458)),
          ("O", (-2.007, 1.439, -0.662))]
  res2 = [("N", (-0.744, 0.404, 0.952)), ("CA", (-0.493, -0.749, 0.119)),
          ("CB", (-1.476, -1.885, 0.331)), ("C", (0.972, -1.193, 0.100)),
          ("O", (1.603, -1.231, 1.129))]
  labels2 = ("B", "A") if swap_second else ("A", "B")
  lines = []
  for resseq, atoms, labels in [(1, res1, ("A", "B")), (2, res2, labels2)]:
    for altloc, shift in zip(labels, (0, 4)):
      for name, (x,y,z) in atoms:
        lines.append(atom_line(name, altloc, "ALA", resseq, (x,y,z+shift),
          name[0]))
  return lines

def exercise_backbone_connectivity():
  h = hierarchy(dipeptide(swap_second=False))
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 0
  assert len(r.unresolved) == 0
  h = hierarchy(dipeptide(swap_second=True))
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 1
  assert len(r.unresolved) == 0
  assert altloc_at(h, 2, (-0.744, 0.404, 0.952)) == "A"
  assert altloc_at(h, 2, (-0.744, 0.404, 4.952)) == "B"

# 1r6j residues 193-196: backbone A/B altlocs, Asp 195 side chain A/C
pdb_str_1r6j_193_196 = """
ATOM      9  N  AALA A 193      -7.535  15.839   9.206  0.82  6.42           N
ATOM     10  N  BALA A 193      -8.074  17.634   9.142  0.18 13.86           N
ATOM     11  CA AALA A 193      -7.488  14.575   8.481  0.82  5.76           C
ATOM     12  CA BALA A 193      -6.837  17.643   8.358  0.18  8.69           C
ATOM     13  C  AALA A 193      -6.094  14.260   7.966  0.82  4.62           C
ATOM     14  C  BALA A 193      -6.660  16.243   7.787  0.18  7.25           C
ATOM     15  O  AALA A 193      -5.774  13.091   7.783  0.82  5.52           O
ATOM     16  O  BALA A 193      -7.265  15.239   8.192  0.18  8.74           O
ATOM     17  CB AALA A 193      -8.454  14.611   7.318  0.82  7.36           C
ATOM     18  CB BALA A 193      -5.717  18.050   9.294  0.18 13.30           C
ATOM     19  N  AMET A 194      -5.267  15.277   7.723  0.82  3.97           N
ATOM     20  N  BMET A 194      -5.739  16.193   6.831  0.18  5.12           N
ATOM     21  CA AMET A 194      -3.954  15.065   7.145  0.82  4.06           C
ATOM     22  CA BMET A 194      -5.129  14.946   6.369  0.18  7.22           C
ATOM     23  C  AMET A 194      -2.841  14.975   8.198  0.82  4.69           C
ATOM     24  C  BMET A 194      -4.143  14.524   7.444  0.18  8.12           C
ATOM     25  O  AMET A 194      -1.659  14.928   7.816  0.82  6.27           O
ATOM     26  O  BMET A 194      -3.120  15.144   7.682  0.18  9.71           O
ATOM     27  CB AMET A 194      -3.639  16.105   6.067  0.82  4.09           C
ATOM     28  CB BMET A 194      -4.453  15.167   5.029  0.18  5.33           C
ATOM     29  CG AMET A 194      -4.671  16.131   4.954  0.82  4.43           C
ATOM     30  CG BMET A 194      -5.295  15.445   3.817  0.18  8.89           C
ATOM     31  SD AMET A 194      -4.741  14.541   4.102  0.82  5.05           S
ATOM     32  SD BMET A 194      -6.566  14.256   3.468  0.18  7.95           S
ATOM     33  CE AMET A 194      -6.433  14.365   3.695  0.82  8.16           C
ATOM     34  CE BMET A 194      -7.923  14.781   4.479  0.18  8.03           C
ATOM     35  N  AASP A 195      -3.175  14.865   9.464  0.82  4.03           N
ATOM     36  N  BASP A 195      -4.530  13.494   8.203  0.18  6.18           N
ATOM     37  CA AASP A 195      -2.145  14.513  10.426  0.82  4.32           C
ATOM     38  CA BASP A 195      -3.700  13.172   9.387  0.18  8.23           C
ATOM     39  C  AASP A 195      -1.602  13.132  10.062  0.82  4.11           C
ATOM     40  C  BASP A 195      -2.677  12.091   9.064  0.18  5.05           C
ATOM     41  O  AASP A 195      -2.384  12.223   9.707  0.82  4.49           O
ATOM     42  O  BASP A 195      -3.017  11.017   8.552  0.18  7.27           O
ATOM     43  CB AASP A 195      -2.668  14.521  11.870  0.82  4.74           C
ATOM     44  CG AASP A 195      -2.873  15.863  12.546  0.41  4.49           C
ATOM     45  CG CASP A 195      -2.549  15.899  12.632  0.41  6.32           C
ATOM     46  OD1AASP A 195      -3.172  16.874  11.866  0.41  5.95           O
ATOM     47  OD1CASP A 195      -2.206  16.861  11.953  0.41  6.00           O
ATOM     48  OD2AASP A 195      -2.747  15.875  13.788  0.41  7.70           O
ATOM     49  OD2CASP A 195      -2.847  15.876  13.835  0.41  7.04           O
ATOM     50  N  APRO A 196      -0.304  12.909  10.205  0.82  4.00           N
ATOM     51  N  BPRO A 196      -1.422  12.324   9.446  0.18  8.89           N
ATOM     52  CA APRO A 196       0.234  11.579   9.921  0.82  3.93           C
ATOM     53  CA BPRO A 196      -0.482  11.195   9.346  0.18  6.96           C
ATOM     54  C  APRO A 196      -0.212  10.593  10.991  0.82  3.42           C
ATOM     55  C  BPRO A 196      -0.779  10.216  10.474  0.18  5.67           C
ATOM     56  O  APRO A 196      -0.685  10.983  12.076  0.82  3.90           O
ATOM     57  O  BPRO A 196      -1.350  10.629  11.482  0.18  7.91           O
ATOM     58  CB APRO A 196       1.765  11.813   9.949  0.82  4.43           C
ATOM     59  CB BPRO A 196       0.863  11.826   9.686  0.18  6.62           C
ATOM     60  CG APRO A 196       1.907  12.950  10.964  0.82  4.79           C
ATOM     61  CG BPRO A 196       0.497  12.961  10.599  0.18  9.15           C
ATOM     62  CD APRO A 196       0.718  13.854  10.683  0.82  4.93           C
ATOM     63  CD BPRO A 196      -0.797  13.506  10.044  0.18 10.44           C
"""

def exercise_real_backbone_altlocs():
  # Consistent labels: nothing to change, no false clashes across peptide
  # bonds (1-3 and 1-4 pairs)
  h = hierarchy(pdb_str_1r6j_193_196.splitlines())
  s = h.as_pdb_string()
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert r.n_sites == 4
  assert len(r.swaps) == 0
  assert len(r.unresolved) == 0
  # Swap A and B of Asp 195: original labels restored, C left alone
  for rg in h.residue_groups():
    if (rg.resseq_as_int() == 195):
      rg.relabel_altlocs({"A":"B", "B":"A"})
  assert h.as_pdb_string() != s
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 1
  assert len(r.unresolved) == 0
  assert h.as_pdb_string() == s

def exercise_adjacent_residues_close_contacts():
  # ARG-GLY from a real structure, B = A shifted by 0.37 A. O of ARG is
  # 2.6 A from C of GLY (five bonds apart): normal geometry, not a clash.
  arg_gly = [
    ("N", "ARG", 434, (-21.288,-36.878,40.267), "N"),
    ("CA", "ARG", 434, (-22.627,-37.426,40.428), "C"),
    ("CB", "ARG", 434, (-23.530,-37.016,39.247), "C"),
    ("CG", "ARG", 434, (-23.748,-35.502,39.151), "C"),
    ("CD", "ARG", 434, (-25.212,-35.145,38.990), "C"),
    ("NE", "ARG", 434, (-25.520,-33.793,39.465), "N"),
    ("CZ", "ARG", 434, (-25.602,-32.710,38.694), "C"),
    ("NH1", "ARG", 434, (-25.391,-32.794,37.386), "N"),
    ("NH2", "ARG", 434, (-25.895,-31.534,39.234), "N"),
    ("C", "ARG", 434, (-22.604,-38.943,40.636), "C"),
    ("O", "ARG", 434, (-23.661,-39.589,40.732), "O"),
    ("N", "GLY", 435, (-21.395,-39.503,40.723), "N"),
    ("CA", "GLY", 435, (-21.231,-40.940,40.825), "C"),
    ("C", "GLY", 435, (-22.324,-41.727,40.118), "C"),
    ("O", "GLY", 435, (-22.346,-41.818,38.889), "O")]
  lines = []
  for resseq in (434, 435):
    for altloc, shift in [("A", (0,0,0)), ("B", (0.3,0.2,0.1))]:
      for name, resname, rs, xyz, e in arg_gly:
        if (rs != resseq): continue
        lines.append(atom_line(name, altloc, resname, resseq,
          [x+s for x, s in zip(xyz, shift)], e))
  h = hierarchy(lines)
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert r.n_sites == 2
  assert len(r.swaps) == 0
  assert len(r.unresolved) == 0

pdb_str_1ee3_his_cd = """\
ATOM   2403  CB AHIS P 303      -0.519  -6.598   3.142  0.50 18.55           C
ATOM   2404  CB BHIS P 303      -0.488  -6.582   3.159  0.50 17.50           C
ATOM   2405  CG AHIS P 303      -1.669  -7.562   3.089  0.50 17.20           C
ATOM   2406  CG BHIS P 303       0.195  -7.565   2.256  0.50 17.42           C
ATOM   2407  ND1AHIS P 303      -2.068  -8.144   1.907  0.50 17.52           N
ATOM   2408  ND1BHIS P 303      -0.068  -8.911   2.333  0.50 17.66           N
ATOM   2409  CD2AHIS P 303      -2.488  -7.976   4.087  0.50 17.55           C
ATOM   2410  CD2BHIS P 303       1.122  -7.345   1.291  0.50 17.48           C
ATOM   2411  CE1AHIS P 303      -3.112  -8.890   2.208  0.50 17.02           C
ATOM   2412  CE1BHIS P 303       0.697  -9.478   1.421  0.50 16.64           C
ATOM   2413  NE2AHIS P 303      -3.409  -8.824   3.512  0.50 17.01           N
ATOM   2414  NE2BHIS P 303       1.436  -8.574   0.763  0.50 17.24           N
HETATM 2431 CD  A CD P 351      -0.316  -9.234   0.858  0.15 15.56          CD
HETATM 2432 CD  B CD P 351      -3.033  -9.002   2.439  0.15 13.13          CD
"""

def exercise_metals():
  # Cd A coordinated by His A (2.3 A), on top of the ring of His B, and
  # vice versa (1ee3): consistent. Mixed pairings put the same ring atoms
  # at bonding distance; metal coordination is not used as a link.
  h = hierarchy(pdb_str_1ee3_his_cd.splitlines())
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 0
  assert len(r.unresolved) == 0
  # Zn A coordinated by SG A, Zn B by SG B, each overlapping the SG of the
  # other conformer. Coordination is not a clash: consistent labels kept,
  # swapped ones restored.
  for swapped in (False, True):
    la, lb = ("B", "A") if swapped else ("A", "B")
    h = hierarchy([
      atom_line("SG", "A", "CYS", 1, (0,0,0), "S"),
      atom_line("SG", "B", "CYS", 1, (3.3,0,0), "S"),
      atom_line("ZN", la, " ZN", 2, (-2.3,0,0), "ZN"),
      atom_line("ZN", lb, " ZN", 2, (1,0,0), "ZN")])
    r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
    assert len(r.swaps) == int(swapped), (swapped, r.swaps)
    assert len(r.unresolved) == 0
    assert altloc_at(h, 2, (-2.3,0,0)) == "A"
  # Consistent Mg coordinated by OP1 of A and B. The other conformer of the
  # nucleotide puts three more O at coordination distance from Mg A: not
  # evidence for relabeling.
  h = hierarchy([
    atom_line("OP1", "A", "  A", 1, (0,0,0), "O"),
    atom_line("OP1", "B", "  A", 1, (10,0,0), "O"),
    atom_line("O5'", "B", "  A", 1, (2.07,2.1,0), "O"),
    atom_line("O3'", "B", "  A", 1, (2.07,-2.1,0), "O"),
    atom_line("O4'", "B", "  A", 1, (2.07,0,2.1), "O"),
    atom_line("MG", "A", " MG", 2, (2.07,0,0), "MG"),
    atom_line("MG", "B", " MG", 2, (12.07,0,0), "MG")])
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 0
  # Chloride is not a metal: at bonding distance it overlaps
  h = hierarchy([
    atom_line("N", "A", "GLY", 1, (0,0,0), "N"),
    atom_line("N", "B", "GLY", 1, (10,0,0), "N"),
    atom_line("CL", "A", " CL", 2, (1.8,0,0), "CL"),
    atom_line("CL", "B", " CL", 2, (20,0,0), "CL")])
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 1
  assert altloc_at(h, 2, (1.8,0,0)) == "B"
  # Chloride with label A only, overlapping conformer B (1r6j): consistent,
  # not a link to B
  h = hierarchy([
    atom_line("N", "A", "GLY", 1, (0,0,0), "N"),
    atom_line("N", "B", "GLY", 1, (10,0,0), "N"),
    atom_line("CL", "A", " CL", 2, (11.8,0,0), "CL")])
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 0

def exercise_links_to_ligands():
  # Covalent link to a ligand (ND2-C1 1.45 A, glycosylation) is a bond, not
  # a clash. Consistent labels are kept, swapped ones restored.
  for swapped in (False, True):
    la, lb = ("B", "A") if swapped else ("A", "B")
    h = hierarchy([
      atom_line("ND2", "A", "ASN", 1, (0,0,0), "N"),
      atom_line("ND2", "B", "ASN", 1, (10,0,0), "N"),
      atom_line("C1", la, "NAG", 2, (1.45,0,0), "C"),
      atom_line("C1", lb, "NAG", 2, (11.45,0,0), "C")])
    r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
    assert len(r.swaps) == int(swapped), (swapped, r.swaps)
    assert len(r.unresolved) == 0
    assert altloc_at(h, 2, (1.45,0,0)) == "A"
  # Asp B and ligand A happen to be at bonding distance: an overlap of
  # alternate conformers, not a link (no matching A-A or B-B bond)
  h = hierarchy([
    atom_line("OD1", "A", "ASP", 1, (0,0,0), "O"),
    atom_line("OD1", "B", "ASP", 1, (10,0,0), "O"),
    atom_line("C5", "A", "ADN", 2, (11.45,0,0), "C"),
    atom_line("C5", "B", "ADN", 2, (30,0,0), "C")])
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 0
  # A ligand atom 0.5 A from a protein atom is an overlap, not a bond
  h = hierarchy([
    atom_line("ND2", "A", "ASN", 1, (0,0,0), "N"),
    atom_line("ND2", "B", "ASN", 1, (10,0,0), "N"),
    atom_line("C1", "A", "NAG", 2, (0.5,0,0), "C"),
    atom_line("C1", "B", "NAG", 2, (20,0,0), "C")])
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 1
  assert altloc_at(h, 2, (0.5,0,0)) == "B"

def exercise_hydrogens():
  # Waters with O and H in each conformer; only the H atoms clash, so the
  # clash is ignored unless include_hydrogens=True
  def w(resseq, altloc, x, o_y, h_y):
    return [atom_line("O", altloc, "HOH", resseq, (x,o_y,0), "O"),
            atom_line("H1", altloc, "HOH", resseq, (x,h_y,0), "H")]
  lines = w(1, "A", 0, 0, 1) + w(1, "B", 8, 0, 1) + \
          w(2, "A", 0, 3.2, 2.2) + w(2, "B", 8, 3.2, 2.2)
  h = hierarchy(lines)
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert len(r.swaps) == 0
  r = h.reconcile_altlocs(crystal_symmetry=None, include_hydrogens=True,
    log=null_out())
  assert len(r.swaps) == 1

def exercise_hd_exchange_labels_kept():
  # Neutron H/D exchange sites are altlocs too (H and D at the same place,
  # hydrogen-only atom groups). Residue 2 has side chain conformers A/B and
  # an exchange site C/D on its blank N: only A and B may be relabeled.
  lines = water(1, [("A", (0,0,0)), ("B", (6,0,0))]) + [
    atom_line("N",  "",  "SER", 2, (3,10,0), "N"),
    atom_line("OG", "A", "SER", 2, (0,2,0),  "O"),
    atom_line("OG", "B", "SER", 2, (6,2,0),  "O"),
    atom_line("H",  "C", "SER", 2, (3,11,0), "H"),
    atom_line("D",  "D", "SER", 2, (3,11,0), "D")]
  for include_hydrogens in [False, True]:
    h = hierarchy(lines)
    r = h.reconcile_altlocs(crystal_symmetry=None,
      include_hydrogens=include_hydrogens, log=null_out())
    assert r.n_sites == 2
    assert len(r.swaps) == 1
    assert r.swaps[0][1] == {"A":"B", "B":"A"}
    assert altloc_at(h, 2, (0,2,0)) == "B"
    for atom in h.atoms():
      if (atom.name.strip() in ["H", "D"]):
        assert atom.parent().altloc == {"H":"C", "D":"D"}[atom.name.strip()]
  # Residues with exchange sites only are not conformers: nothing to do
  h = hierarchy([
    atom_line("H", "A", "SER", 1, (0,0,0), "H"),
    atom_line("D", "B", "SER", 1, (0,0,0), "D"),
    atom_line("H", "A", "SER", 2, (0,1.2,0), "H"),
    atom_line("D", "B", "SER", 2, (0,1.2,0), "D")])
  r = h.reconcile_altlocs(crystal_symmetry=None, include_hydrogens=True,
    log=null_out())
  assert r.n_sites == 0
  assert len(r.swaps) == 0

def exercise_greedy_search():
  # Large groups of sites are searched greedily: same results as the
  # exhaustive search on small cases
  from iotbx.pdb import altloc_reconciliation
  h_backbone = hierarchy(pdb_str_1r6j_193_196.splitlines())
  for rg in h_backbone.residue_groups():
    if (rg.resseq_as_int() == 195):
      rg.relabel_altlocs({"A":"B", "B":"A"})
  # Row of waters 2 A apart, each clashing with its neighbors unless labels
  # alternate sides; waters 3 and 4 both swapped: relabeling either alone
  # removes one clash and makes another
  water_row = []
  for i in range(6):
    sides = [(0,2*i,0), (6,2*i,0)]
    if (i % 2 == 1 or i in (2, 3)) and not (i % 2 == 1 and i in (2, 3)):
      sides.reverse()
    water_row += water(i+1, [("A", sides[0]), ("B", sides[1])])
  # two adjacent residues swapped: no single relabeling helps
  h_backbone_2 = hierarchy(pdb_str_1r6j_193_196.splitlines())
  for rg in h_backbone_2.residue_groups():
    if (rg.resseq_as_int() in (194, 195)):
      rg.relabel_altlocs({"A":"B", "B":"A"})
  cases = [
    two_inconsistent_waters(),
    water(1, [("A", (0,0,0)), ("B", (6,0,0))]) +
    water(2, [("A", (0,2,0)), ("B", (6,2,0))]) +
    water(3, [("A", (0,4,0)), ("B", (6,4,0))]),
    water(1, [("A", (0,0,0)), ("B", (3.6,0,0)), ("C", (0,10,0))]) +
    water(2, [("A", (1.8,1,0)), ("B", (5.6,0,0)), ("C", (0,20,0))]),
    water(1, [("A", (0,0,0)),   ("B", (20,0,0))]) +
    water(2, [("A", (2,0,0)),   ("B", (22,0,0))]) +
    water(3, [("A", (1,1.7,0)), ("B", (21,1.7,0))]),
    water(1, [("A", (0,0,0)), ("B", (10,0,0))]) +
    water(2, [("A", (2,0,0))]) +
    water(3, [("A", (12,0,0))]),
    h_backbone.as_pdb_string().splitlines(),
    h_backbone_2.as_pdb_string().splitlines(),
    water_row]
  for lines in cases:
    results = []
    for max_combinations in [10000, 1]:
      h = hierarchy(lines)
      r = altloc_reconciliation.run(hierarchy=h, crystal_symmetry=None,
        max_combinations=max_combinations)
      results.append((r, h))
    (r_exhaustive, h_exhaustive), (r_greedy, h_greedy) = results
    assert r_greedy.unresolved == r_exhaustive.unresolved
    # greedy may pick another, equally consistent labeling (e.g. swap the
    # neighbors instead): nothing left to change for the exhaustive search
    r = altloc_reconciliation.run(hierarchy=h_greedy, crystal_symmetry=None)
    assert len(r.swaps) == 0, r.swaps
    if (len(r_exhaustive.unresolved) == 0):
      assert len(r_greedy.swaps) == len(r_exhaustive.swaps)

def exercise_no_altlocs():
  h = hierarchy(water(1, [("", (0,0,0))]) + water(2, [("", (0,2,0))]))
  r = h.reconcile_altlocs(crystal_symmetry=None, log=null_out())
  assert r.n_sites == 0
  assert len(r.swaps) == 0

def run():
  exercise_relabel_altlocs()
  exercise_two_waters_swapped()
  exercise_check_only()
  exercise_distant_sites_untouched()
  exercise_fewest_changes()
  exercise_three_labels()
  exercise_unresolvable_flagged()
  exercise_unresolvable_kept_local()
  exercise_marginal_frustration_ignored()
  exercise_marginal_clash_ignored()
  exercise_insignificant_gain_ignored()
  exercise_symmetry_contact()
  exercise_nonsense_symmetry()
  exercise_clash_through_other_conformer()
  exercise_geometry_clash_not_flagged()
  exercise_partial_labels_flagged()
  exercise_backbone_connectivity()
  exercise_real_backbone_altlocs()
  exercise_adjacent_residues_close_contacts()
  exercise_metals()
  exercise_links_to_ligands()
  exercise_hydrogens()
  exercise_hd_exchange_labels_kept()
  exercise_greedy_search()
  exercise_no_altlocs()
  print("OK")

if (__name__ == "__main__"):
  run()
