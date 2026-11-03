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

if __name__ == "__main__":
  test_side_chain_split()
  test_one_conformer()
  test_conformers_disagree()
  test_extra_sigma_bonds()
  test_run_without_restraints()
  print("OK")
