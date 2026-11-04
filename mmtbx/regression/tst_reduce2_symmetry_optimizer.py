"""reduce2 qmready=True on a hydroxyl H at a crystal contact.

Fixture: DG A 12 and DG B 22 of 1BNA with the 1BNA cell. O3' of DG A 12 (a
3' end) accepts an H-bond from H22 of DG B 22 across a crystal contact
(symmetry mate -x+3/2,-y+1,z-1/2). qmready places HO3' on it, and turned
about C3'-O3' against the asymmetric unit alone it lands 0.56 A from that
H22, although other orientations are 1.79 A from everything.

Two things keep that from the Sorry "Bad clash(es) prevents qmready=True"
of the symmetry-aware clash check (cutoff 1.0 A):
- include_symmetric_copies=True scores the Movers against the symmetry mates
  as a fixed environment, so the optimizer itself turns HO3' clear of H22.
- With the default (no copies) the qmready regularization pushes the clash
  apart until the violation check is satisfied, so HO3' ends just over the
  cutoff and the run succeeds.
Both are asserted; the second only at the cutoff, since the regularization
stops at the first violation-free step.
"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
from iotbx.cli_parser import run_program
from libtbx.utils import null_out
from mmtbx.programs import reduce2
from scitbx import matrix

cryst1 = "CRYST1   24.870   40.390   66.200  90.00  90.00  90.00 P 21 21 21\n"

# 1BNA DG A 12 and DG B 22 (waters removed)
pdb_two_dg = """\
ATOM    222  P    DG A  12      19.331  25.774  -7.925  1.00 55.98           P
ATOM    223  OP1  DG A  12      20.704  25.976  -7.408  1.00 45.83           O
ATOM    224  OP2  DG A  12      18.763  26.851  -8.758  1.00 44.26           O
ATOM    225  O5'  DG A  12      19.302  24.412  -8.763  1.00 62.63           O
ATOM    226  C5'  DG A  12      20.109  23.284  -8.359  1.00 69.50           C
ATOM    227  C4'  DG A  12      19.748  22.167  -9.299  1.00 39.92           C
ATOM    228  O4'  DG A  12      18.350  21.969  -9.139  1.00 32.00           O
ATOM    229  C3'  DG A  12      19.921  22.404 -10.815  1.00 50.39           C
ATOM    230  O3'  DG A  12      20.985  21.635 -11.401  1.00 64.13           O
ATOM    231  C2'  DG A  12      18.535  22.062 -11.381  1.00 36.18           C
ATOM    232  C1'  DG A  12      17.965  21.200 -10.269  1.00 24.79           C
ATOM    233  N9   DG A  12      16.493  21.220 -10.265  1.00 28.44           N
ATOM    234  C8   DG A  12      15.663  22.289 -10.478  1.00 31.85           C
ATOM    235  N7   DG A  12      14.368  21.958 -10.390  1.00 38.26           N
ATOM    236  C5   DG A  12      14.388  20.640 -10.102  1.00 28.99           C
ATOM    237  C6   DG A  12      13.301  19.742  -9.856  1.00 42.63           C
ATOM    238  O6   DG A  12      12.091  19.967  -9.857  1.00 49.17           O
ATOM    239  N1   DG A  12      13.750  18.466  -9.625  1.00 40.15           N
ATOM    240  C2   DG A  12      15.042  18.043  -9.605  1.00 33.42           C
ATOM    241  N2   DG A  12      15.259  16.717  -9.406  1.00 40.53           N
ATOM    242  N3   DG A  12      16.061  18.885  -9.792  1.00 37.34           N
ATOM    243  C4   DG A  12      15.660  20.156 -10.027  1.00 31.14           C
TER
ATOM    424  P    DG B  22       6.509  21.324  20.099  1.00 56.50           P
ATOM    425  OP1  DG B  22       5.387  20.397  20.396  1.00 50.81           O
ATOM    426  OP2  DG B  22       6.235  22.774  20.306  1.00 53.84           O
ATOM    427  O5'  DG B  22       7.767  20.924  20.993  1.00 66.30           O
ATOM    428  C5'  DG B  22       8.216  19.559  21.073  1.00 73.42           C
ATOM    429  C4'  DG B  22       9.422  19.557  21.977  1.00 42.96           C
ATOM    430  O4'  DG B  22      10.493  20.260  21.319  1.00 52.87           O
ATOM    431  C3'  DG B  22       9.267  20.267  23.325  1.00 38.51           C
ATOM    432  O3'  DG B  22      10.088  19.657  24.293  1.00 60.28           O
ATOM    433  C2'  DG B  22       9.751  21.670  22.990  1.00 22.00           C
ATOM    434  C1'  DG B  22      10.988  21.226  22.256  1.00 24.85           C
ATOM    435  N9   DG B  22      11.599  22.357  21.543  1.00 25.91           N
ATOM    436  C8   DG B  22      11.037  23.545  21.159  1.00 23.91           C
ATOM    437  N7   DG B  22      11.921  24.362  20.566  1.00 39.18           N
ATOM    438  C5   DG B  22      13.072  23.653  20.580  1.00 25.66           C
ATOM    439  C6   DG B  22      14.370  24.003  20.102  1.00 28.34           C
ATOM    440  O6   DG B  22      14.747  25.057  19.585  1.00 31.85           O
ATOM    441  N1   DG B  22      15.268  22.983  20.308  1.00 25.22           N
ATOM    442  C2   DG B  22      15.023  21.776  20.891  1.00 11.07           C
ATOM    443  N2   DG B  22      16.066  20.914  21.038  1.00 25.92           N
ATOM    444  N3   DG B  22      13.815  21.452  21.350  1.00 19.05           N
ATOM    445  C4   DG B  22      12.902  22.429  21.151  1.00 23.69           C
TER
END
"""

CLASH_CUTOFF = 1.0  # the qmready clash check

def run(text, prefix, extra_args=()):
  fn = prefix + ".pdb"
  with open(fn, "w") as f:
    f.write(text)
  out = prefix + "_out.pdb"
  if os.path.exists(out): os.remove(out)
  return run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn, "qmready=True", "overwrite=True", "output.filename=" + out] +
      list(extra_args))

def atom_index(hierarchy, chain_id, resseq):
  """{name: i_seq} of one residue"""
  result = {}
  for atom in hierarchy.atoms():
    rg = atom.parent().parent()
    if rg.resseq.strip() == resseq and rg.parent().id.strip() == chain_id:
      result[atom.name.strip()] = atom.i_seq
  return result

def min_symmetry_distance(hierarchy, crystal_symmetry, i_seq, exclude):
  """Shortest distance from atom i_seq to any atom not in exclude, over all
  symmetry mates (4.0 if none that close), and that atom's id with the
  symmetry operation"""
  xrs = hierarchy.extract_xray_structure(crystal_symmetry=crystal_symmetry)
  table = xrs.pair_asu_table(distance_cutoff=4.0).extract_pair_sym_table()
  cell = crystal_symmetry.unit_cell()
  sites_frac = xrs.sites_frac()
  atoms = hierarchy.atoms()
  best = (4.0, None)
  for i, j, invert in [(i_seq, j, False) for j in range(len(table))] + [
                      (i, i_seq, True) for i in range(len(table))]:
    other = j if not invert else i
    if other == i_seq or other in exclude: continue
    for op in table[i].get(j, []):
      d = cell.distance(sites_frac[i], op * sites_frac[j])
      if d < best[0]:
        best = (d, "%s %s" % (atoms[other].id_str(), op))
  return best

def best_rotation(hierarchy, crystal_symmetry, idx):
  """Turn HO3' of the residue idx about C3'-O3' in 10 degree steps and
  return the largest achievable min_symmetry_distance to atoms outside the
  residue (and the angle); the hierarchy is left as it was"""
  atoms = hierarchy.atoms()
  i_h, i_o, i_c = idx["HO3'"], idx["O3'"], idx["C3'"]
  o = matrix.col(atoms[i_o].xyz)
  axis = (o - matrix.col(atoms[i_c].xyz)).normalize()
  v0 = matrix.col(atoms[i_h].xyz) - o
  exclude = set(idx.values())
  best = (-1, None)
  for deg in range(0, 360, 10):
    r = axis.axis_and_angle_as_r3_rotation_matrix(deg, deg=True)
    atoms[i_h].set_xyz(o + r * v0)
    d = min_symmetry_distance(hierarchy, crystal_symmetry, i_h, exclude)[0]
    if d > best[0]: best = (d, deg)
  atoms[i_h].set_xyz(o + v0)
  return best

def test_hydroxyl_at_crystal_contact():
  prefix = "tst_reduce2_symmetry_optimizer"
  cs = iotbx.pdb.input(source_info=None, lines=cryst1).crystal_symmetry()
  # Without CRYST1 qmready places HO3' on the 3' end against the asymmetric
  # unit alone: it sits on H22 of the symmetry mate of DG B 22, although
  # another orientation about C3'-O3' is clear of everything.
  r = run(pdb_two_dg, prefix + "_p1")
  h = r.model.get_hierarchy()
  idx = atom_index(h, "A", "12")
  assert "HO3'" in idx, sorted(idx)
  d_best, deg = best_rotation(h, cs, idx)
  assert d_best >= CLASH_CUTOFF, (d_best, deg)
  d_placed, who = min_symmetry_distance(h, cs, idx["HO3'"], set(idx.values()))
  print("HO3' DG A 12 as placed without the cell: %.2f A from %s; best "
        "orientation %.2f A at %d degrees" % (d_placed, who, d_best, deg))
  assert d_placed < CLASH_CUTOFF, (d_placed, who)
  # With the cell and the symmetry mates as the optimizer's environment,
  # HO3' is turned clear of them by a margin.
  r = run(cryst1 + pdb_two_dg, prefix + "_copies",
    extra_args=["include_symmetric_copies=True"])
  h = r.model.get_hierarchy()
  idx = atom_index(h, "A", "12")
  d, who = min_symmetry_distance(h, cs, idx["HO3'"], set(idx.values()))
  print("with include_symmetric_copies=True: %.2f A from %s" % (d, who))
  assert d >= 1.5, (d, who)
  # With the defaults the regularization clears the clash check.
  r = run(cryst1 + pdb_two_dg, prefix)
  h = r.model.get_hierarchy()
  idx = atom_index(h, "A", "12")
  d, who = min_symmetry_distance(h, cs, idx["HO3'"], set(idx.values()))
  print("with the defaults: %.2f A from %s" % (d, who))
  assert d >= CLASH_CUTOFF, (d, who)

if __name__ == "__main__":
  test_hydroxyl_at_crystal_contact()
  print("OK")
