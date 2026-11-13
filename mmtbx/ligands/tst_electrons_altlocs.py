"""mmtbx.ligands.electrons on a model with alternate conformations.

electron_distribution counts the bonds of the restraints manager, where a
shared atom (altloc blank) is bonded to its neighbour in every conformer:
CA of a residue whose side chain is split is bonded to CB A and CB B. Counted
on the whole model those bonds make the shared atoms charged and the total
charge wrong (the fragment below, -1, came out as -4 from reduce2). The
distribution is therefore built one conformer at a time (altloc blank + X),
mapped back to the i_seqs of the whole model, and the total charge is the
charge the conformers agree on.
"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out, Sorry
from mmtbx.ligands import electrons

# reduce2 qmready=True on ASP A 5 - GLY A 6 with the side chain (and HA) of
# the ASP in altlocs A and B. NH3+ and COO-, ASP -1: charge -1.
pdb_str = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM     20  N   ASP A   5       9.580   7.432   8.699  1.00 37.74           N
ATOM     21  CA  ASP A   5      10.692   8.369   8.630  1.00 38.33           C
ATOM     22  C   ASP A   5      10.828   9.231   9.884  1.00 37.17           C
ATOM     23  O   ASP A   5      11.938   9.511  10.322  1.00 36.50           O
ATOM     24  H1  ASP A   5       9.676   6.724   7.944  1.00 37.74           H
ATOM     25  H2  ASP A   5       9.585   6.956   9.624  1.00 37.74           H
ATOM     26  H3  ASP A   5       8.685   7.947   8.579  1.00 37.74           H
ATOM     27  CB AASP A   5      10.625   9.197   7.354  0.50 39.43           C
ATOM     28  CG AASP A   5      11.281   8.467   6.194  0.50 43.72           C
ATOM     29  OD1AASP A   5      10.549   7.831   5.406  0.50 49.92           O
ATOM     30  OD2AASP A   5      12.531   8.382   6.067  0.50 51.10           O
ATOM     31  HA AASP A   5      11.619   7.797   8.584  0.50 38.33           H
ATOM     32  HB2AASP A   5       9.582   9.385   7.098  0.50 39.43           H
ATOM     33  HB3AASP A   5      11.147  10.141   7.507  0.50 39.43           H
ATOM     34  CB BASP A   5      10.925   9.197   7.354  0.50 39.43           C
ATOM     35  CG BASP A   5      11.581   8.467   6.194  0.50 43.72           C
ATOM     36  OD1BASP A   5      10.849   7.831   5.406  0.50 49.92           O
ATOM     37  OD2BASP A   5      12.831   8.382   6.067  0.50 51.10           O
ATOM     38  HA BASP A   5      11.529   7.674   8.560  0.50 38.33           H
ATOM     39  HB2BASP A   5       9.959   9.558   7.001  0.50 39.43           H
ATOM     40  HB3BASP A   5      11.568  10.040   7.608  0.50 39.43           H
ATOM     41  N   GLY A   6       9.691   9.552  10.506  1.00 36.16           N
ATOM     42  CA  GLY A   6       9.668  10.262  11.762  1.00 35.61           C
ATOM     43  C   GLY A   6      10.374   9.397  12.850  1.00 34.61           C
ATOM     44  O   GLY A   6      11.228   9.895  13.566  1.00 33.43           O
ATOM     45  OXT GLY A   6      10.070   8.212  12.987  1.00 33.43           O
ATOM     46  H   GLY A   6       8.763   9.325  10.150  1.00 36.16           H
ATOM     47  HA2 GLY A   6      10.191  11.213  11.664  1.00 35.61           H
ATOM     48  HA3 GLY A   6       8.639  10.452  12.066  1.00 35.61           H
END
"""

def lines(drop=(), keep_altloc=None, altloc_of=None):
  """pdb_str without the atoms named in drop ((name, altloc) pairs); with
  keep_altloc only altloc blank and keep_altloc, the label blanked; with
  altloc_of ({name: altloc}) those blank atoms moved into that altloc"""
  out = []
  for l in pdb_str.splitlines():
    if l.startswith("ATOM"):
      name, altloc = l[12:16].strip(), l[16]
      if (name, altloc) in drop:
        continue
      if keep_altloc is not None:
        if altloc not in (" ", keep_altloc):
          continue
        l = l[:16] + " " + l[17:]
      if altloc_of and name in altloc_of and altloc == " ":
        l = l[:16] + altloc_of[name] + l[17:]
    out.append(l)
  return "\n".join(out) + "\n"

def get_model(text):
  inp = iotbx.pdb.input(lines=text.splitlines(), source_info=None)
  model = mmtbx.model.manager(model_input=inp, log=null_out())
  model.process(make_restraints=True)
  return model

def distribution(text):
  return electrons.run(model=get_model(text), return_formal_charges=True)

def i_seq(model, resseq, name, altloc=" "):
  for atom in model.get_hierarchy().atoms():
    if (atom.parent().parent().resseq.strip() == str(resseq) and
        atom.name.strip() == name and
        atom.parent().altloc.strip() == altloc.strip()):
      return atom.i_seq
  raise AssertionError((resseq, name, altloc))

def formal_charges(ed, altloc=" "):
  """{(resseq, name): electrons} of the atoms of altloc blank + altloc"""
  rc = {}
  for atom in ed.hierarchy.atoms():
    ag = atom.parent()
    if ag.altloc.strip() in ("", altloc.strip()):
      rc[(ag.parent().resseq, atom.name.strip())] = ed[atom.i_seq]
  return rc

def test_side_chain_split():
  """-1; every conformer gives the formal charges of the same conformer
  without altlocs; the atoms of the whole model are the keys"""
  ed = distribution(pdb_str)
  assert ed.get_total_charge() == -1, ed.get_total_charge()
  assert ed.get_conformer_charges() == {"A": -1, "B": -1}, \
    ed.get_conformer_charges()
  atoms = ed.hierarchy.atoms()
  assert sorted(k for k in ed if type(k) != tuple) == list(range(len(atoms)))
  for altloc in "AB":
    ref = distribution(lines(keep_altloc=altloc))
    assert ref.get_total_charge() == -1
    assert ref.get_conformer_charges() == {}, ref.get_conformer_charges()
    assert formal_charges(ed, altloc) == formal_charges(ref), altloc
  # NH3+ (-1 electron); one O of each carboxylate +1
  fc = formal_charges(ed, "A")
  assert fc[("   5", "N")] == -1, fc
  assert sorted(v for k, v in fc.items() if v) == [-1, 1, 1], fc
  # the shared atom bonded to both conformers, CA bonded to CB A, CB B,
  # HA A, HA B, is not charged, nor are the split atoms
  assert fc[("   5", "CA")] == 0, fc
  for altloc in "AB":
    for name in ("CB", "CG", "HA"):
      assert fc[("   5", name)] == 0, (altloc, name, fc)

def test_one_conformer():
  """alternative_location_id=X: the distribution of that conformer only"""
  model = get_model(pdb_str)
  for altloc in "AB":
    ed = electrons.electron_distribution(
      model.get_hierarchy(), model.get_restraints_manager().geometry,
      alternative_location_id=altloc)
    assert ed.get_total_charge() == -1, ed.get_total_charge()
    assert ed.get_conformer_charges() == {altloc: -1}
    present = sorted(k for k in ed if type(k) != tuple)
    expected = sorted(a.i_seq for a in ed.hierarchy.atoms()
                      if a.parent().altloc.strip() in ("", altloc))
    assert present == expected, (altloc, len(present), len(expected))

def test_conformers_disagree():
  """H3 of ASP A 5 in altloc A only: NH3+ in A (-1), NH2 in B (-2);
  no single total charge"""
  ed = distribution(lines(altloc_of={"H3": "A"}))
  assert ed.get_conformer_charges() == {"A": -1, "B": -2}, \
    ed.get_conformer_charges()
  try:
    ed.get_total_charge()
  except Sorry as e:
    assert "A: -1" in str(e) and "B: -2" in str(e), str(e)
  else:
    raise AssertionError("no Sorry")

def test_extra_sigma_bonds():
  """the capping path: HA2 of GLY A 6 and HB3 of both ASP conformers
  removed and replaced with one extra sigma bond on CA, CB A and CB B"""
  text = lines(drop=[("HA2", " "), ("HB3", "A"), ("HB3", "B")])
  model = get_model(text)
  extra = {i_seq(model, 6, "CA"): 1, i_seq(model, 5, "CB", "A"): 1,
           i_seq(model, 5, "CB", "B"): 1}
  ed = electrons.run(model=model, return_formal_charges=True,
                     extra_sigma_bonds=extra)
  assert ed.get_total_charge() == -1, ed.get_total_charge()
  for i in extra:
    assert ed[i] == 0, (i, ed[i])
  # without the extra bonds the CB of each conformer is open-shell
  ed = electrons.run(model=model, return_formal_charges=True)
  for altloc in "AB":
    i = i_seq(model, 5, "CB", altloc)
    assert ed[i] == 1, (altloc, ed[i])

def test_run_without_restraints():
  """the legacy entry point of qr.finalise: hierarchy and symmetry"""
  inp = iotbx.pdb.input(lines=pdb_str.splitlines(), source_info=None)
  charge = electrons.run(pdb_hierarchy=inp.construct_hierarchy(),
                         crystal_symmetry=inp.crystal_symmetry())
  assert charge == -1, charge


def with_altloc(text, resseq, altloc):
  """pdb_str with the blank atoms of residue resseq moved into altloc"""
  out = []
  for l in text.splitlines():
    if l.startswith("ATOM") and l[22:26].strip() == str(resseq) and l[16] == " ":
      l = l[:16] + altloc + l[17:]
    out.append(l)
  return "\n".join(out) + "\n"

def test_inherited_conformer():
  """GLY A 6 in altloc C only, ASP A 5 in A and B: a residue without atoms
  of a conformer belongs to it through its first own altloc (as
  pdb_interpretation inherits it along the chain), so every conformer is
  the whole fragment and the charges agree; without that rule conformer C
  held the ASP as its blank atoms only and came out -4"""
  ed = distribution(with_altloc(pdb_str, 6, "C"))
  assert ed.get_conformer_charges() == {"A": -1, "B": -1, "C": -1}, \
    ed.get_conformer_charges()
  assert ed.get_total_charge() == -1
  sizes = dict((a, d.hierarchy.atoms_size()) for a, d in ed.conformers.items())
  assert len(set(sizes.values())) == 1, sizes
  for altloc, d in ed.conformers.items():
    # the ASP of conformer C is conformer A's
    expected = "A" if altloc == "C" else altloc
    names = set(a.name.strip() for a in d.hierarchy.atoms())
    assert "CB" in names and "OXT" in names, (altloc, names)
    for atom in d.hierarchy.atoms():
      if atom.name.strip() == "CB":
        assert abs(atom.xyz[0] - (10.625 if expected == "A" else 10.925)) < 1e-3

# 1ejg THR A 21 - PRO/SER A 22 - GLU A 23 as deposited (with H): residue 22 is
# PRO in altloc A and SER in B and C with one blank N, which iotbx reads as
# two residue groups (PRO: blank N and A; SER: B and C).
two_groups_str = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM    386  N   THR A  21       4.056   9.404  -4.216  1.00  1.67
ATOM    387  CA  THR A  21       3.664  10.576  -3.461  1.00  1.43
ATOM    388  C   THR A  21       4.915  11.347  -3.001  1.00  1.58
ATOM    389  O   THR A  21       5.844  10.729  -2.478  1.00  2.25
ATOM    390  CB  THR A  21       2.932  10.113  -2.178  1.00  1.52
ATOM    391  OG1 THR A  21       1.964   9.124  -2.552  1.00  1.89
ATOM    392  CG2 THR A  21       2.280  11.274  -1.459  1.00  2.12
ATOM    393  H   THR A  21       3.676   8.507  -3.952  1.00  2.00
ATOM    394  HA  THR A  21       3.014  11.225  -4.064  1.00  1.71
ATOM    395  HB  THR A  21       3.666   9.650  -1.503  1.00  1.83
ATOM    396  HG1 THR A  21       2.092   8.329  -2.016  1.00  2.85
ATOM    397 HG21 THR A  21       1.606  11.800  -2.148  1.00  3.17
ATOM    398 HG22 THR A  21       3.056  11.968  -1.105  1.00  3.17
ATOM    399 HG23 THR A  21       1.706  10.898  -0.600  1.00  3.17
ATOM    400  N   PRO A  22       4.915  12.683  -3.102  1.00  1.83
ATOM    401  CA APRO A  22       6.042  13.429  -2.601  0.57  1.82
ATOM    402  C  APRO A  22       6.387  13.122  -1.160  0.57  1.66
ATOM    403  O  APRO A  22       5.480  13.006  -0.345  0.57  2.09
ATOM    404  CB APRO A  22       5.655  14.896  -2.744  0.57  2.86
ATOM    405  CG APRO A  22       4.661  14.854  -4.058  0.57  2.66
ATOM    406  CD APRO A  22       3.957  13.505  -3.910  0.57  2.27
ATOM    407  HA APRO A  22       6.915  13.244  -3.241  0.57  1.93
ATOM    408  HB2APRO A  22       5.112  15.266  -1.863  0.57  3.34
ATOM    409  HB3APRO A  22       6.516  15.539  -2.944  0.57  3.34
ATOM    410  HG2APRO A  22       3.919  15.683  -4.043  0.57  3.20
ATOM    411  HG3APRO A  22       5.214  14.903  -5.002  0.57  3.20
ATOM    412  HD2APRO A  22       2.984  13.625  -3.389  0.57  2.69
ATOM    413  HD3APRO A  22       3.752  13.060  -4.903  0.57  2.69
ATOM    414  CA BSER A  22       6.034  13.399  -2.687  0.21  1.55
ATOM    415  CA CSER A  22       6.112  13.653  -2.656  0.22  2.31
ATOM    416  C  BSER A  22       6.367  13.062  -1.223  0.21  2.92
ATOM    417  C  CSER A  22       6.354  13.275  -1.187  0.22  1.92
ATOM    418  O  BSER A  22       5.412  13.050  -0.345  0.21  1.87
ATOM    419  O  CSER A  22       5.636  12.705  -0.270  0.22  1.77
ATOM    420  CB BSER A  22       5.409  14.835  -2.876  0.21  2.60
ATOM    421  CB CSER A  22       5.605  15.097  -2.687  0.22  3.56
ATOM    422  OG BSER A  22       4.760  15.243  -1.635  0.21  2.11
ATOM    423  OG CSER A  22       6.750  15.771  -2.280  0.22  6.38
ATOM    424  HA BSER A  22       6.897  13.240  -3.349  0.21  4.26
ATOM    425  HA CSER A  22       7.013  13.518  -3.271  0.22  4.26
ATOM    426  N   GLU A  23       7.651  13.190  -0.860  1.00  1.86
ATOM    427  CA  GLU A  23       8.108  13.002   0.513  1.00  1.81
ATOM    428  C   GLU A  23       7.394  13.968   1.459  1.00  1.66
ATOM    429  O   GLU A  23       7.064  13.585   2.589  1.00  1.97
ATOM    430  CB  GLU A  23       9.620  13.168   0.681  1.00  2.12
ATOM    431  CG AGLU A  23      10.459  12.148  -0.048  0.58  2.73
ATOM    432  CG BGLU A  23      10.414  11.969   0.088  0.42  1.79
ATOM    433  CD AGLU A  23      11.942  12.208   0.242  0.58  4.60
ATOM    434  CD BGLU A  23      11.911  12.048   0.452  0.42  2.33
ATOM    435  OE1AGLU A  23      12.418  13.330   0.541  0.58  7.57
ATOM    436  OE1BGLU A  23      12.430  13.113   0.795  0.42  2.70
ATOM    437  OE2AGLU A  23      12.629  11.211   0.128  0.58  7.42
ATOM    438  OE2BGLU A  23      12.514  10.945   0.316  0.42  4.43
ATOM    439  H   GLU A  23       8.361  13.404  -1.574  1.00  2.23
ATOM    440  HA  GLU A  23       7.842  11.977   0.808  1.00  2.17
ATOM    441  HB2 GLU A  23       9.902  14.169   0.326  1.00  2.55
ATOM    442  HB3 GLU A  23       9.860  13.116   1.752  1.00  2.55
ATOM    443  HG2AGLU A  23      10.094  11.146   0.214  0.58  3.27
ATOM    444  HG2BGLU A  23       9.996  11.029   0.473  0.42  2.13
ATOM    445  HG3AGLU A  23      10.310  12.283  -1.128  0.58  3.27
ATOM    446  HG3BGLU A  23      10.306  11.966  -1.006  0.42  2.13
END
"""

def test_two_residue_groups():
  """the residue groups of one residue are taken together: conformer A is
  not given SER B (the first own altloc of the second group), whose CA
  would be a fourth bond on the shared N; all conformers agree"""
  ed = distribution(two_groups_str)
  charges = ed.get_conformer_charges()
  assert len(set(charges.values())) == 1, charges
  assert ed.get_total_charge() == charges["A"]
  assert ed.get_conformer_charge_differences() == [], \
    ed.get_conformer_charge_differences()
  for altloc, d in ed.conformers.items():
    resnames = set(a.parent().resname for a in d.hierarchy.atoms()
                   if a.parent().parent().resseq.strip() == "22")
    assert resnames == set(["PRO" if altloc == "A" else "SER", "PRO"]), \
      (altloc, resnames)
    # the shared N: C, CA A and CD A in conformer A. In B and C the H of
    # the SER is in the SER group and pdb_interpretation does not bond it
    # to the N of the other group (a lone H): N -1 and H +1 cancel
    if altloc != "A": continue
    for atom in d.hierarchy.atoms():
      if atom.name.strip() == "N" and atom.parent().resname == "PRO":
        assert d[atom.i_seq] == 0, (altloc, d[atom.i_seq])

# 1pq7 CYS A 41 (A, B) and CYS A 57 (blank), with the H reduce2 places: the
# disulfide SG A 41 - SG 57 takes HG A 57, HG B 57 stays (an H-only atom
# group B on the blank SG). Two residues cut from the chain: the same
# truncation charges in both conformers.
cys_str = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM     25  N  ACYS A  41       5.508   7.360   1.752  0.59  3.97           N
ATOM     26  CA ACYS A  41       4.721   6.164   1.517  0.60  3.70           C
ATOM     27  C  ACYS A  41       3.333   6.308   2.146  0.60  3.80           C
ATOM     28  O  ACYS A  41       2.979   7.330   2.757  0.60  3.74           O
ATOM     29  CB ACYS A  41       5.431   4.898   2.010  0.61  3.50           C
ATOM     30  SG ACYS A  41       6.731   4.335   0.827  0.60  3.82           S
ATOM     31  H  ACYS A  41       5.397   7.692   2.538  0.59  3.97           H
ATOM     32  HA ACYS A  41       4.609   6.046   0.561  0.60  3.70           H
ATOM     33  HB2ACYS A  41       5.855   5.081   2.863  0.61  3.50           H
ATOM     34  HB3ACYS A  41       4.780   4.185   2.109  0.61  3.50           H
ATOM     35  N  BCYS A  41       5.374   7.503   1.660  0.40  4.85           N
ATOM     36  CA BCYS A  41       4.595   6.308   1.311  0.39  5.11           C
ATOM     37  C  BCYS A  41       3.340   6.305   2.209  0.40  4.25           C
ATOM     38  O  BCYS A  41       3.124   7.213   3.049  0.40  3.20           O
ATOM     39  CB BCYS A  41       5.483   5.140   1.450  0.40  7.12           C
ATOM     40  SG BCYS A  41       5.836   4.613   3.064  0.39 10.10           S
ATOM     41  H  BCYS A  41       5.164   7.892   2.398  0.40  4.85           H
ATOM     42  HA BCYS A  41       4.272   6.271   0.397  0.39  5.11           H
ATOM     43  HB2BCYS A  41       5.068   4.392   0.992  0.40  7.12           H
ATOM     44  HB3BCYS A  41       6.331   5.356   1.031  0.40  7.12           H
ATOM     45  HG BCYS A  41       6.022   3.427   3.056  0.39 10.10           H
ATOM     71  N   CYS A  57       9.445   0.516   1.315  1.00  4.44           N
ATOM     72  CA  CYS A  57       9.313   1.815   1.946  1.00  4.72           C
ATOM     73  C   CYS A  57      10.241   1.988   3.148  1.00  4.55           C
ATOM     74  O   CYS A  57      10.553   3.132   3.509  1.00  5.19           O
ATOM     75  CB  CYS A  57       7.875   2.020   2.439  1.00  5.19           C
ATOM     76  SG  CYS A  57       6.674   2.294   1.132  1.00  6.27           S
ATOM     77  H   CYS A  57       8.749   0.015   1.378  1.00  4.44           H
ATOM     78  HA  CYS A  57       9.546   2.481   1.280  1.00  4.72           H
ATOM     79  HB2 CYS A  57       7.600   1.229   2.929  1.00  5.19           H
ATOM     80  HB3 CYS A  57       7.856   2.795   3.022  1.00  5.19           H
ATOM     81  HG BCYS A  57       7.054   3.193   0.434  0.39  6.27           H
END
"""

def test_h_only_group():
  """An H-only atom group is not a conformer of its residue: conformer A of
  CYS 57 is the blank atoms (SG bonded to SG A 41), without HG B."""
  ed = distribution(cys_str)
  charges = ed.get_conformer_charges()
  assert charges["A"] == charges["B"], charges
  assert ed.get_conformer_charge_differences() == [], \
    ed.get_conformer_charge_differences()
  n = dict((a, ed.conformers[a].atoms.size()) for a in ed.conformers)
  # B has HG B 57 and HG B 41, A neither; the heavy atoms are the same count
  assert n["B"] == n["A"] + 2, n

def test_differences_named():
  """the Sorry of get_total_charge() names the residues that differ"""
  ed = distribution(lines(altloc_of={"H3": "A"}))
  # the residue's own charge: NH3+ and COO- in A, NH2 and COO- in B
  assert ed.get_conformer_charge_differences() == \
    [("A 5", {"A": 0, "B": -1})], ed.get_conformer_charge_differences()
  try:
    ed.get_total_charge()
  except Sorry as e:
    assert "Residues whose charge differs between the conformers: A 5 " \
      "(A: 0, B: -1)" in str(e), str(e)
  else:
    raise AssertionError("no Sorry")
if __name__ == "__main__":
  test_side_chain_split()
  test_one_conformer()
  test_conformers_disagree()
  test_extra_sigma_bonds()
  test_run_without_restraints()
  test_inherited_conformer()
  test_two_residue_groups()
  test_h_only_group()
  test_differences_named()
  print("OK")
