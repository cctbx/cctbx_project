from __future__ import absolute_import, division, print_function
"""
backrub_move: the flanking C-alphas stay put, the moved block is rigid, the
CA-CB direction turns by the angle; with the peptide counter-rotations the
bond lengths are kept and both carbonyl oxygens end closer to where they
were. Residues given as atom groups or as conformer residues.
"""
import math
import iotbx.pdb
from scitbx import matrix
from mmtbx.refinement.real_space.fit_residue import backrub_move

pdb_str = """\
ATOM    440  N   LEU A  59      -7.510 122.749  38.503  1.00 89.42           N
ATOM    441  CA  LEU A  59      -8.480 121.839  37.894  1.00 90.13           C
ATOM    442  C   LEU A  59      -9.914 122.306  38.123  1.00 90.59           C
ATOM    443  O   LEU A  59     -10.778 122.134  37.254  1.00 88.28           O
ATOM    444  CB  LEU A  59      -8.304 120.428  38.434  1.00 82.68           C
ATOM    445  CG  LEU A  59      -9.142 119.474  37.595  1.00 86.20           C
ATOM    446  CD1 LEU A  59      -8.562 119.366  36.180  1.00 93.60           C
ATOM    447  CD2 LEU A  59      -9.250 118.118  38.258  1.00 93.32           C
ATOM    448  N   PHE A  60     -10.187 122.893  39.292  1.00 89.57           N
ATOM    449  CA  PHE A  60     -11.493 123.501  39.534  1.00 89.43           C
ATOM    450  C   PHE A  60     -11.793 124.577  38.497  1.00 88.98           C
ATOM    451  O   PHE A  60     -12.907 124.648  37.970  1.00 87.16           O
ATOM    452  CB  PHE A  60     -11.538 124.079  40.950  1.00 86.28           C
ATOM    453  CG  PHE A  60     -12.899 124.556  41.382  1.00 85.30           C
ATOM    454  CD1 PHE A  60     -13.947 123.661  41.533  1.00 88.93           C
ATOM    455  CD2 PHE A  60     -13.120 125.892  41.665  1.00 81.10           C
ATOM    456  CE1 PHE A  60     -15.193 124.092  41.937  1.00 81.71           C
ATOM    457  CE2 PHE A  60     -14.358 126.329  42.074  1.00 80.56           C
ATOM    458  CZ  PHE A  60     -15.397 125.427  42.209  1.00 83.74           C
ATOM    459  N   GLU A  61     -10.801 125.413  38.183  1.00 87.78           N
ATOM    460  CA  GLU A  61     -10.971 126.431  37.152  1.00 89.55           C
ATOM    461  C   GLU A  61     -11.309 125.807  35.805  1.00 93.53           C
ATOM    462  O   GLU A  61     -12.214 126.275  35.104  1.00 93.94           O
ATOM    463  CB  GLU A  61      -9.703 127.270  37.029  1.00 89.15           C
ATOM    464  CG  GLU A  61      -9.646 128.501  37.893  1.00 81.08           C
ATOM    465  CD  GLU A  61      -8.344 129.250  37.693  1.00 92.28           C
ATOM    466  OE1 GLU A  61      -7.298 128.739  38.156  1.00 87.99           O
ATOM    467  OE2 GLU A  61      -8.362 130.324  37.045  1.00 95.81           O
"""

BONDS = [(0, "C", 1, "N"), (1, "N", 1, "CA"), (1, "CA", 1, "C"), (1, "C", 1, "O"),
         (1, "CA", 1, "CB"), (1, "C", 2, "N"), (0, "C", 0, "O")]

def residues(use_atom_groups):
  h = iotbx.pdb.input(source_info=None, lines=pdb_str).construct_hierarchy()
  rgs = list(h.only_model().chains()[0].residue_groups())
  if(use_atom_groups): return h, [rg.atom_groups()[0] for rg in rgs]
  return h, [rg.conformers()[0].only_residue() for rg in rgs]

def xyz(res, name):
  for a in res.atoms():
    if(a.name.strip() == name): return matrix.col(a.xyz)

def dist(rs, b):
  return abs(xyz(rs[b[0]], b[1]) - xyz(rs[b[2]], b[3]))

def exercise(use_atom_groups):
  for move_oxygens in (False, True):
    h0, r0 = residues(use_atom_groups)
    h, r = residues(use_atom_groups)
    backrub_move(prev_res=r[0], cur_res=r[1], next_res=r[2], angle=10.,
      move_oxygens=move_oxygens, accept_worse_rama=True)
    for k in (0, 2):
      assert abs(xyz(r[k], "CA") - xyz(r0[k], "CA")) < 1.e-9
    for b in BONDS:
      assert abs(dist(r, b) - dist(r0, b)) < 1.e-6, (move_oxygens, b,
        dist(r, b), dist(r0, b))
    v0 = xyz(r0[1], "CB") - xyz(r0[1], "CA")
    v1 = xyz(r[1], "CB") - xyz(r[1], "CA")
    turn = math.degrees(math.acos(min(1., v0.dot(v1) / abs(v0) / abs(v1))))
    assert 8. < turn < 10.5, turn
    shift_o = [abs(xyz(r[k], "O") - xyz(r0[k], "O")) for k in (0, 1)]
    if(not move_oxygens): rigid = shift_o
    else:
      for s, s_rigid in zip(shift_o, rigid):
        assert s < s_rigid, (shift_o, rigid)

if(__name__ == "__main__"):
  exercise(use_atom_groups=True)
  exercise(use_atom_groups=False)
  print("OK")
