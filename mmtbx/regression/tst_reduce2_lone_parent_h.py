"""reduce2 places the H of a parent without heavy neighbours (NH4+): riding
H cannot parameterize them, so they were dropped ("riding H refused") and the
ion was written as a bare N, counted as N3- by the charge calculation. Now
they are placed as free H (mmtbx.hydrogens.lone_parent_h): ideal geometry,
oriented by H-bonds and clearance, reported under their own heading, and the
charge of 2PVB-like models is right (NH4+ counts +1)."""
from __future__ import absolute_import, division, print_function
import math
import iotbx.pdb
from scitbx import matrix
from scitbx.array_family import flex
from mmtbx.hydrogens import lone_parent_h
from mmtbx.regression.tst_reduce2_capping import run, pdb_lys_asp_gly

# ASP A 5 - GLY A 6 of pdb_lys_asp_gly and an ammonium ion 2.83 A from OD1
# of ASP 5 as read (10.549 7.831 5.406), away from CG
nh4 = """\
HETATM   20  N   NH4 A 200       8.900   6.400   3.600  1.00 15.00           N
"""
pdb_str = "".join(l + "\n" for l in pdb_lys_asp_gly.splitlines()
  if not (l.startswith("ATOM") and l[17:20] == "LYS")).replace("END\n",
  nh4 + "END\n")

def atoms_of(hierarchy, resseq):
  return dict((a.name.strip(), a) for a in hierarchy.atoms()
              if a.parent().parent().resseq.strip() == str(resseq))

def check_nh4(hierarchy, nh_length):
  a = atoms_of(hierarchy, 200)
  n = a["N"]
  hs = [a[k] for k in sorted(a) if k != "N"]
  assert [h.name.strip() for h in hs] == ["HN1", "HN2", "HN3", "HN4"], \
    sorted(a)
  for h in hs:
    assert abs(h.distance(n) - nh_length) < 0.01, (h.name, h.distance(n))
    assert h.occ == n.occ, (h.occ, n.occ)
  for i in range(4):
    for j in range(i + 1, 4):
      u = matrix.col(hs[i].xyz) - matrix.col(n.xyz)
      v = matrix.col(hs[j].xyz) - matrix.col(n.xyz)
      assert abs(u.angle(v, deg=True) - 109.47) < 1.0, u.angle(v, deg=True)
  # one H donates to the nearer carboxylate O of ASP 5 (OD1 as read; reduce2
  # may rename OD1/OD2 to the IUPAC convention), straight
  asp = atoms_of(hierarchy, 5)
  acceptor = min((asp["OD1"], asp["OD2"]), key=lambda o: o.distance(n))
  assert abs(acceptor.distance(n) - 2.83) < 0.05, acceptor.distance(n)
  best = min(hs, key=lambda h: h.distance(acceptor))
  assert best.distance(acceptor) < 2.2, best.distance(acceptor)
  u = matrix.col(n.xyz) - matrix.col(best.xyz)
  v = matrix.col(acceptor.xyz) - matrix.col(best.xyz)
  assert u.angle(v, deg=True) > 150, u.angle(v, deg=True)
  # clear of everything else
  for h in hs:
    for b in hierarchy.atoms():
      if b.parent().parent().resseq.strip() == "200": continue
      d = h.distance(b)
      if b.element_is_hydrogen(): assert d > 1.9, (h.name, b.id_str(), d)
      else: assert d > 1.5 or b is acceptor, (h.name, b.id_str(), d)

def test_add():
  prefix = "tst_reduce2_lone_parent_1"
  fn_out = prefix + "_out.pdb"
  run(pdb_str, prefix, [], fn_out=fn_out)
  h = iotbx.pdb.input(fn_out).construct_hierarchy()
  check_nh4(h, 0.89)
  with open(prefix + "_out.txt") as f:
    description = f.read()
  assert "riding H refused" not in description, description
  assert "H placed as free atoms" in description and \
    "NH4 A 200: HN1 HN2 HN3 HN4 on N" in description, description

def test_qmready_charge():
  prefix = "tst_reduce2_lone_parent_2"
  fn_out = prefix + "_out.pdb"
  r = run(pdb_str, prefix, ["qmready=True"], fn_out=fn_out)
  h = iotbx.pdb.input(fn_out).construct_hierarchy()
  check_nh4(h, 1.04)
  # ASP-GLY -1 (NH3+, COO-, ASP), NH4+ +1
  assert r.total_charge == 0, r.total_charge
  with open(prefix + "_out.txt") as f:
    assert "qmready: H not placed" not in f.read()

def test_isolated_group_geometry():
  # no environment: a tetrahedron at the bond length, deterministic
  sites = lone_parent_h.place((1.0, 2.0, 3.0), [1.0] * 4, 109.47,
    environment=lone_parent_h.group_args(xyz=flex.vec3_double(),
      clearance=flex.double(), acceptor=[]))
  assert len(sites) == 4
  p = matrix.col((1.0, 2.0, 3.0))
  for s in sites:
    assert abs((matrix.col(s) - p).length() - 1.0) < 1e-6
  again = lone_parent_h.place((1.0, 2.0, 3.0), [1.0] * 4, 109.47,
    environment=lone_parent_h.group_args(xyz=flex.vec3_double(),
      clearance=flex.double(), acceptor=[]))
  assert sites == again
  # two and three H: the requested angle
  for n, angle in ((2, 104.5), (3, 107.0)):
    d = lone_parent_h.ideal_directions(n, angle)
    for i in range(n):
      for j in range(i + 1, n):
        assert abs(d[i].angle(d[j], deg=True) - angle) < 1e-6

if __name__ == "__main__":
  test_add()
  test_qmready_charge()
  test_isolated_group_geometry()
  print("OK")
