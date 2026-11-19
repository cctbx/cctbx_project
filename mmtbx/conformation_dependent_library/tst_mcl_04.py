"""
Symmetry contacts in the Fe-S hook, get_sulfur_iron_cluster_coordination(): a
cluster on a symmetry element or with a ligand across a contact gets no
restraints (b576e99c42); any other contact is skipped; only that cluster counts.
Fragments: 6O0B, 7TPZ, 5L7L as deposited; 9IYU with a water placed by hand.
"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
from cctbx.geometry_restraints import linking_class
from libtbx.utils import null_out

origin_id = linking_class.linking_class().get_origin_id('metal coordination')

# SF4 with its four Cys; water A 503 of a mate 3.34 A from S1
pdb_6o0b = """\
CRYST1   57.071   93.045   60.796  90.00  98.88  90.00 P 1 21 1
ATOM    647  N   CYS A  97      -1.477  -0.038  62.988  1.00 20.86           N
ATOM    648  CA  CYS A  97      -2.583  -0.388  63.861  1.00 20.94           C
ATOM    649  C   CYS A  97      -3.457  -1.416  63.156  1.00 20.01           C
ATOM    650  O   CYS A  97      -3.756  -1.282  61.966  1.00 21.38           O
ATOM    651  CB  CYS A  97      -3.378   0.864  64.232  1.00 19.93           C
ATOM    652  SG  CYS A  97      -4.929   0.477  65.101  1.00 21.65           S
ATOM    901  N   CYS A 132      -8.098  -4.942  71.864  1.00 20.14           N
ATOM    902  CA  CYS A 132      -8.079  -3.486  71.777  1.00 20.57           C
ATOM    903  C   CYS A 132      -9.047  -3.009  70.691  1.00 19.17           C
ATOM    904  O   CYS A 132      -9.600  -3.806  69.930  1.00 20.15           O
ATOM    905  CB  CYS A 132      -6.649  -3.018  71.472  1.00 21.59           C
ATOM    906  SG  CYS A 132      -6.031  -3.660  69.893  1.00 20.01           S
ATOM   2689  N   CYS B  97      -6.857   1.711  73.173  1.00 22.69           N
ATOM   2690  CA  CYS B  97      -5.545   1.416  72.604  1.00 20.94           C
ATOM   2691  C   CYS B  97      -4.552   1.206  73.738  1.00 25.86           C
ATOM   2692  O   CYS B  97      -4.447   2.041  74.644  1.00 22.71           O
ATOM   2693  CB  CYS B  97      -5.095   2.547  71.677  1.00 20.11           C
ATOM   2694  SG  CYS B  97      -3.413   2.394  71.037  1.00 20.10           S
ATOM   2943  N   CYS B 132       2.675  -3.360  67.345  1.00 18.86           N
ATOM   2944  CA  CYS B 132       1.975  -2.276  66.680  1.00 21.04           C
ATOM   2945  C   CYS B 132       2.481  -0.931  67.202  1.00 23.44           C
ATOM   2946  O   CYS B 132       3.265  -0.856  68.153  1.00 21.39           O
ATOM   2947  CB  CYS B 132       0.455  -2.404  66.910  1.00 19.79           C
ATOM   2948  SG  CYS B 132       0.005  -2.375  68.660  1.00 20.74           S
HETATM 4131 FE1  SF4 A 301      -4.606  -1.932  69.299  1.00 19.63          FE
HETATM 4132 FE2  SF4 A 301      -3.485   0.520  69.745  1.00 19.68          FE
HETATM 4133 FE3  SF4 A 301      -4.129  -0.328  67.143  1.00 19.43          FE
HETATM 4134 FE4  SF4 A 301      -2.082  -1.480  68.550  1.00 18.69          FE
HETATM 4135  S1  SF4 A 301      -2.134   0.707  67.817  1.00 20.37           S
HETATM 4136  S2  SF4 A 301      -3.659  -2.574  67.265  1.00 19.75           S
HETATM 4137  S3  SF4 A 301      -2.870  -1.427  70.760  1.00 21.25           S
HETATM 4138  S4  SF4 A 301      -5.647   0.118  68.880  1.00 19.51           S
HETATM 4241  O   HOH A 503     -17.662 -43.194  54.107  1.00 25.70           O
END
"""

# half an SF4 on a 2-fold, its own mate within 3.5 A; two Cys
pdb_7tpz = """\
CRYST1   45.522   73.711   74.665  90.00  90.00  90.00 P 2 21 21
ATOM    701  N   CYS A  97      15.132 -34.734  -5.151  1.00 28.38           N
ATOM    702  CA  CYS A  97      14.272 -34.627  -3.981  1.00 23.79           C
ATOM    703  C   CYS A  97      13.141 -33.644  -4.283  1.00 22.58           C
ATOM    704  O   CYS A  97      13.395 -32.506  -4.691  1.00 22.66           O
ATOM    705  CB  CYS A  97      15.104 -34.171  -2.783  1.00 24.81           C
ATOM    706  SG  CYS A  97      14.190 -33.647  -1.342  1.00 23.40           S
ATOM    972  N   CYS A 132       9.993 -34.023   5.254  1.00 20.44           N
ATOM    973  CA  CYS A 132      11.401 -34.317   5.037  1.00 20.25           C
ATOM    974  C   CYS A 132      11.902 -33.411   3.931  1.00 23.71           C
ATOM    975  O   CYS A 132      11.117 -32.802   3.204  1.00 22.45           O
ATOM    976  CB  CYS A 132      11.608 -35.790   4.685  1.00 21.45           C
ATOM    977  SG  CYS A 132      10.506 -36.394   3.417  1.00 24.48           S
HETATM 2081 FE1  SF4 A 301      13.270 -35.566  -0.481  1.00 21.13          FE
HETATM 2082 FE2  SF4 A 301      11.473 -36.485   1.323  1.00 21.83          FE
HETATM 2083  S3  SF4 A 301      10.992 -35.121  -0.449  1.00 21.89           S
HETATM 2084  S4  SF4 A 301      13.747 -36.248   1.648  1.00 20.66           S
END
"""

# half an FES (occupancy 0.51) on a 2-fold; the Cys of its other Fe in a mate
pdb_5l7l = """\
CRYST1   73.760  161.410   93.220  90.00  90.00  90.00 C 2 2 21
ATOM    123  N   CYS A  27     -31.832   2.209   6.355  1.00 92.73           N
ATOM    124  CA  CYS A  27     -32.516   2.743   5.180  1.00 81.73           C
ATOM    125  C   CYS A  27     -34.006   2.883   5.476  1.00 79.53           C
ATOM    126  O   CYS A  27     -34.674   1.899   5.782  1.00 88.59           O
ATOM    127  CB  CYS A  27     -32.288   1.838   3.962  1.00 79.07           C
ATOM    128  SG  CYS A  27     -33.278   2.236   2.492  1.00 95.22           S
ATOM    149  N   CYS A  30     -36.036   0.573   3.268  1.00 83.19           N
ATOM    150  CA  CYS A  30     -35.875  -0.873   3.442  1.00 80.76           C
ATOM    151  C   CYS A  30     -37.035  -1.498   4.218  1.00 80.95           C
ATOM    152  O   CYS A  30     -37.604  -0.862   5.104  1.00 82.48           O
ATOM    153  CB  CYS A  30     -34.551  -1.182   4.146  1.00 73.07           C
ATOM    154  SG  CYS A  30     -33.131  -1.278   3.041  1.00102.60           S
HETATM 3079 FE1  FES A 502     -33.217  -0.252  -1.324  0.51 72.47          FE
HETATM 3080 FE2  FES A 502     -33.217   0.252   1.324  0.51 72.47          FE
HETATM 3081  S1  FES A 502     -31.481   0.000   0.000  0.51 78.51           S
HETATM 3082  S2  FES A 502     -34.958   0.000   0.000  0.51 65.57           S
END
"""

# FES on a 2-fold, its own mate within 3.5 A; SF4 with four Cys, water of a mate
# (placed by hand) 3.3 A from S1 and nothing else
pdb_9iyu = """\
CRYST1  112.638   76.281   72.084  90.00 112.36  90.00 C 1 2 1
ATOM    291  N   CYS A  39       4.473 -23.455   4.413  1.00 28.85           N
ATOM    292  CA  CYS A  39       4.317 -23.538   2.966  1.00 27.70           C
ATOM    293  C   CYS A  39       5.593 -23.148   2.236  1.00 29.15           C
ATOM    294  O   CYS A  39       6.130 -22.060   2.452  1.00 28.69           O
ATOM    295  CB  CYS A  39       3.185 -22.620   2.519  1.00 26.48           C
ATOM    296  SG  CYS A  39       3.013 -22.768   0.738  1.00 25.71           S
ATOM    344  N   CYS A  47       4.480 -17.385  -4.743  1.00 21.68           N
ATOM    345  CA  CYS A  47       3.023 -17.387  -4.586  1.00 21.10           C
ATOM    346  C   CYS A  47       2.488 -15.975  -4.780  1.00 23.90           C
ATOM    347  O   CYS A  47       3.126 -15.001  -4.362  1.00 26.09           O
ATOM    348  CB  CYS A  47       2.583 -17.912  -3.215  1.00 26.72           C
ATOM    349  SG  CYS A  47       0.780 -17.763  -3.048  1.00 25.27           S
ATOM    350  N   CYS A  48       1.332 -15.861  -5.434  1.00 23.28           N
ATOM    351  CA  CYS A  48       0.753 -14.549  -5.709  1.00 21.19           C
ATOM    352  C   CYS A  48      -0.763 -14.625  -5.584  1.00 25.37           C
ATOM    353  O   CYS A  48      -1.408 -15.551  -6.112  1.00 21.29           O
ATOM    354  CB  CYS A  48       1.150 -14.030  -7.090  1.00 21.57           C
ATOM    355  SG  CYS A  48       0.398 -12.401  -7.572  1.00 23.04           S
ATOM    377  N   CYS A  51      -4.850 -11.038  -7.469  1.00 24.50           N
ATOM    378  CA  CYS A  51      -6.138 -10.910  -8.126  1.00 25.02           C
ATOM    379  C   CYS A  51      -6.952 -12.170  -7.838  1.00 23.43           C
ATOM    380  O   CYS A  51      -6.431 -13.170  -7.331  1.00 23.64           O
ATOM    381  CB  CYS A  51      -5.973 -10.652  -9.630  1.00 23.77           C
ATOM    382  SG  CYS A  51      -5.800 -12.105 -10.719  1.00 24.85           S
ATOM    411  N   CYS A  56      -2.087 -18.081  -7.323  1.00 24.86           N
ATOM    412  CA  CYS A  56      -1.137 -18.699  -8.224  1.00 24.98           C
ATOM    413  C   CYS A  56       0.123 -19.075  -7.463  1.00 26.42           C
ATOM    414  O   CYS A  56       0.553 -18.387  -6.524  1.00 24.68           O
ATOM    415  CB  CYS A  56      -0.767 -17.765  -9.368  1.00 25.19           C
ATOM    416  SG  CYS A  56      -2.143 -17.471 -10.453  1.00 26.03           S
ATOM    518  N   CYS A  70      -1.866 -15.095 -15.241  1.00 26.60           N
ATOM    519  CA  CYS A  70      -1.828 -13.951 -16.158  1.00 25.22           C
ATOM    520  C   CYS A  70      -0.582 -13.943 -17.024  1.00 28.81           C
ATOM    521  O   CYS A  70      -0.448 -13.054 -17.869  1.00 25.43           O
ATOM    522  CB  CYS A  70      -1.863 -12.624 -15.403  1.00 25.15           C
ATOM    523  SG  CYS A  70      -0.493 -12.452 -14.252  1.00 24.05           S
HETATM 4673 FE1  SF4 A 701      -1.331 -12.999 -12.062  1.00 24.01          FE
HETATM 4674 FE2  SF4 A 701      -1.090 -12.953  -9.338  1.00 21.27          FE
HETATM 4675 FE3  SF4 A 701      -2.175 -15.066 -10.562  1.00 24.58          FE
HETATM 4676 FE4  SF4 A 701      -3.549 -12.766 -10.428  1.00 25.27          FE
HETATM 4677  S1  SF4 A 701      -2.973 -14.035  -8.666  1.00 22.99           S
HETATM 4678  S2  SF4 A 701      -3.304 -14.151 -12.250  1.00 25.44           S
HETATM 4679  S3  SF4 A 701      -1.829 -11.241 -10.588  1.00 22.18           S
HETATM 4680  S4  SF4 A 701      -0.017 -14.311 -10.853  1.00 23.26           S
HETATM 4681 FE1  FES A 702       0.464 -19.257  -1.220  1.00 28.03          FE
HETATM 4682 FE2  FES A 702       1.306 -21.188   0.386  1.00 27.46          FE
HETATM 4683  S1  FES A 702       0.654 -21.421  -1.762  1.00 27.70           S
HETATM 4684  S2  FES A 702       1.884 -19.028   0.551  1.00 26.09           S
HETATM 9001  O   HOH A 901       5.459 -15.937   9.710  1.00 20.00           O
END
"""

def proxies(pdb_str):
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.process(make_restraints=True)
  atoms = model.get_hierarchy().atoms()
  grm = model.get_restraints_manager().geometry
  bond_proxies = grm.pair_proxies(sites_cart=atoms.extract_xyz()).bond_proxies
  bonds = []
  for p in bond_proxies.simple:
    if p.origin_id == origin_id:
      bonds.append((atoms[p.i_seqs[0]].id_str(), atoms[p.i_seqs[1]].id_str(),
        p.distance_ideal, p.weight))
  for p in bond_proxies.asu:
    assert p.origin_id != origin_id
  angles = []
  for p in grm.angle_proxies:
    if p.origin_id == origin_id:
      angles.append((tuple(atoms[i].id_str() for i in p.i_seqs), p.angle_ideal,
        p.weight))
  return sorted(bonds), sorted(angles)

def without(pdb_str, resnames):
  return "\n".join(l for l in pdb_str.split("\n") if l[17:20] not in resnames)

def test_other_contact_is_skipped():
  # rule: a contact that is not the cluster's own mate or a ligand is skipped
  bonds, angles = proxies(pdb_6o0b)
  assert len(bonds) == 4, bonds
  assert len(angles) == 12, angles
  assert proxies(without(pdb_6o0b, ["HOH"])) == (bonds, angles)

def test_cluster_on_a_symmetry_element():
  # rule: a cluster against its own mate gets no restraints
  assert proxies(pdb_7tpz) == ([], [])

def test_ligand_across_a_contact():
  # rule: a cluster with a ligand across a contact gets no restraints
  assert proxies(pdb_5l7l) == ([], [])

def test_only_that_cluster():
  # rule: the first two apply per cluster: FES none, SF4 its own
  bonds, angles = proxies(pdb_9iyu)
  assert len(bonds) == 4, bonds
  assert len(angles) == 12, angles
  for a1, a2, ideal, weight in bonds:
    assert "FES" not in a1 + a2 and "SF4" in a1 + a2
  assert proxies(without(pdb_9iyu, ["HOH", "FES"])) == (bonds, angles)

if (__name__ == "__main__"):
  test_other_contact_is_skipped()
  test_cluster_on_a_symmetry_element()
  test_ligand_across_a_contact()
  test_only_that_cluster()
  print("OK")
