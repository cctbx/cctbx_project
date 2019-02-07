from __future__ import division, print_function
import time
import mmtbx.model
import iotbx.pdb
import mmtbx.hydrogens
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal


def run():
  exercise1()
  exercise2()
  exercise3()
  exercise4()
  exercise5()
  exercise6()

# ------------------------------------------------------------------------------

def get_model(pdb_str):
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input = pdb_inp,
                              log         = null_out())
  model_with_h = mmtbx.hydrogens.add(model = model)
  return model_with_h

# ------------------------------------------------------------------------------
# make sure that H of incomplete peptide unit (N-terminal) is not placed
# test will need to be adapted once H3 can be placed at N-terminal
# ------------------------------------------------------------------------------
def exercise1():
  model_with_h = get_model(pdb_str = pdb_str_1)
  number_h = model_with_h.get_hd_selection().count(True)
  # total number of H atoms in Tyr: 9, but peptide H atom cannot be parameterized
  assert(number_h == 8)

# ------------------------------------------------------------------------------
# if one heavy atom is missing, it can impact parameterization and thus
# placement of SEVERAL H atoms.
# ------------------------------------------------------------------------------
def exercise2():
  model_with_h = get_model(pdb_str = pdb_str_2)
  number_h = model_with_h.get_hd_selection().count(True)
  geometry_restraints = model_with_h.get_restraints_manager().geometry
  e_sites = geometry_restraints.energies_sites(
          sites_cart        = model_with_h.get_sites_cart(),
          compute_gradients = False)
  t = e_sites.target
  # The following atoms cannot be parameterized: H, HD1, HE1
  # number_h = 9 - 3 = 6
  assert(number_h == 6)
  # if the target is higher, indicative that H atom is at wrong place
  #assert approx_equal(t,145, 10)
  assert(t < 180)

# residual as of 2/6/2019
#    target: 145.107
#      bond_residual_sum (n=16): 0.00813544
#      nonbonded_residual_sum (n=50): 139.635
#      angle_residual_sum (n=24): 0.304979
#      dihedral_residual_sum (n=5): 5.1542
#      chirality_residual_sum (n=1): 0.000916602
#      planarity_residual_sum (n=1): 0.0031465
#      parallelity_residual_sum (n=0): 0

# with old version of connectivity (less strict for existence of dihedral angles),
# the HE1 atom was placed incorrectly, which yielded following residual
#    target: 343.706
#      bond_residual_sum (n=18): 0.00813544
#      nonbonded_residual_sum (n=76): 286.773
#      angle_residual_sum (n=26): 0.304979
#      dihedral_residual_sum (n=5): 5.1542
#      chirality_residual_sum (n=1): 0.000916602
#      planarity_residual_sum (n=1): 51.4654
#      parallelity_residual_sum (n=0): 0
# new version of connectivity ignores HE1 --> it is not parameterized

# ------------------------------------------------------------------------------
# Don't put H on disulfides
# ------------------------------------------------------------------------------
def exercise3():
  model_with_h = get_model(pdb_str = pdb_str_3)
  hd_sel = model_with_h.get_hd_selection()
  number_h = hd_sel.count(True)
  h_atoms = model_with_h.get_hierarchy().select(hd_sel).atoms()
  h_names = list(h_atoms.extract_name())

  assert(number_h == 6)
  assert('HG' not in h_names)

# ------------------------------------------------------------------------------
# Put H on lone Cys
# ------------------------------------------------------------------------------
def exercise4():
  model_with_h = get_model(pdb_str = pdb_str_4)
  hd_sel = model_with_h.get_hd_selection()
  number_h = hd_sel.count(True)
  h_atoms = model_with_h.get_hierarchy().select(hd_sel).atoms()
  h_names = list(h_atoms.extract_name())

  assert(number_h == 4)
  assert('HG' in h_names)

# ------------------------------------------------------------------------------
# N atom is missing, make sure peptide H is not placed
# ------------------------------------------------------------------------------
def exercise5():
  model_with_h = get_model(pdb_str = pdb_str_5)
  hd_sel = model_with_h.get_hd_selection()
  number_h = hd_sel.count(True)
  h_atoms = model_with_h.get_hierarchy().select(hd_sel).atoms()
  h_names = list(h_atoms.extract_name())

  assert(number_h == 8)
  assert('H' not in h_names)

# ------------------------------------------------------------------------------
# Don't put H on Cys coordinating a metal
# ------------------------------------------------------------------------------
def exercise6():
  model_with_h = get_model(pdb_str = pdb_str_6)
  hd_sel = model_with_h.get_hd_selection()
  number_h = hd_sel.count(True)
  h_atoms = model_with_h.get_hierarchy().select(hd_sel).atoms()
  h_names = list(h_atoms.extract_name())

  print(h_names)
  assert(number_h == 3)
  assert('HG' not in h_names)

# ------------------------------------------------------------------------------

pdb_str_1 = """
CRYST1   14.074   14.385   15.410  90.00  90.00  90.00 P 1
SCALE1      0.071053  0.000000  0.000000        0.00000
SCALE2      0.000000  0.069517  0.000000        0.00000
SCALE3      0.000000  0.000000  0.064893        0.00000
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
"""

pdb_str_2 = """
CRYST1   17.955   13.272   13.095  90.00  90.00  90.00 P 1
ATOM   1942  N   TYR A 139      10.241   7.920   5.000  1.00 10.00           N
ATOM   1943  CA  TYR A 139      10.853   7.555   6.271  1.00 10.00           C
ATOM   1944  C   TYR A 139      12.362   7.771   6.227  1.00 10.00           C
ATOM   1945  O   TYR A 139      12.955   8.272   7.181  1.00 10.00           O
ATOM   1946  CB  TYR A 139      10.540   6.098   6.617  1.00 10.00           C
ATOM   1947  CG  TYR A 139       9.063   5.805   6.749  1.00 10.00           C
ATOM   1949  CD2 TYR A 139       8.414   5.943   7.969  1.00 10.00           C
ATOM   1950  CE1 TYR A 139       6.966   5.122   5.770  1.00 10.00           C
ATOM   1951  CE2 TYR A 139       7.064   5.676   8.095  1.00 10.00           C
ATOM   1952  CZ  TYR A 139       6.345   5.266   6.993  1.00 10.00           C
ATOM   1953  OH  TYR A 139       5.000   5.000   7.113  1.00 10.00           O
TER
"""

pdb_str_3 = """
CRYST1   14.197   15.507   16.075  90.00  90.00  90.00 P 1
ATOM      1  N   CYS A 128      14.558 -30.378 -19.478  1.00 11.73           N
ATOM      2  CA  CYS A 128      13.499 -29.511 -19.901  1.00 11.51           C
ATOM      3  C   CYS A 128      12.454 -30.316 -20.640  1.00 11.99           C
ATOM      4  O   CYS A 128      12.765 -31.329 -21.265  1.00 16.86           O
ATOM      5  CB  CYS A 128      14.059 -28.357 -20.691  1.00 12.94           C
ATOM      6  SG  CYS A 128      15.213 -27.350 -19.760  1.00 12.36           S
ATOM      7  N   CYS A 232      13.105 -25.822 -15.371  1.00  4.75           N
ATOM      8  CA  CYS A 232      14.298 -26.646 -15.573  1.00  5.50           C
ATOM      9  C   CYS A 232      15.596 -25.838 -15.551  1.00  5.70           C
ATOM     10  O   CYS A 232      16.651 -26.346 -15.190  1.00  8.25           O
ATOM     11  CB  CYS A 232      14.178 -27.497 -16.824  1.00  6.69           C
ATOM     12  SG  CYS A 232      14.060 -26.477 -18.360  1.00  7.70           S
TER
"""

pdb_str_4 = """
CRYST1   14.197   15.507   16.075  90.00  90.00  90.00 P 1
ATOM      1  N   CYS A 128      14.558 -30.378 -19.478  1.00 11.73           N
ATOM      2  CA  CYS A 128      13.499 -29.511 -19.901  1.00 11.51           C
ATOM      3  C   CYS A 128      12.454 -30.316 -20.640  1.00 11.99           C
ATOM      4  O   CYS A 128      12.765 -31.329 -21.265  1.00 16.86           O
ATOM      5  CB  CYS A 128      14.059 -28.357 -20.691  1.00 12.94           C
ATOM      6  SG  CYS A 128      15.213 -27.350 -19.760  1.00 12.36           S
TER
"""

pdb_str_5 = """
CRYST1   15.000   15.000   15.000  90.00 90.00  90.00 P 1           1
ATOM      1  CB  PHE H   1       7.767   5.853   7.671  1.00 15.00           C
ATOM      2  CG  PHE H   1       6.935   5.032   8.622  1.00 15.00           C
ATOM      3  CD1 PHE H   1       5.918   4.176   8.140  1.00 15.00           C
ATOM      4  CD2 PHE H   1       7.161   5.107  10.012  1.00 15.00           C
ATOM      5  CE1 PHE H   1       5.126   3.395   9.038  1.00 15.00           C
ATOM      6  CE2 PHE H   1       6.382   4.336  10.930  1.00 15.00           C
ATOM      7  CZ  PHE H   1       5.360   3.476  10.439  1.00 15.00           C
ATOM      8  C   PHE H   1       7.956   7.811   6.133  1.00 15.00           C
ATOM      9  O   PHE H   1       8.506   7.237   5.169  1.00 15.00           O
ATOM     15  CA  PHE H   1       7.000   7.000   7.000  1.00 15.00           C
END
"""

pdb_str_6 = """
CRYST1  108.910  108.910  108.910  90.00  90.00  90.00 I 4 3 2      48
ATOM    366  N   CYS A  48     -40.561 -16.057   6.830  1.00 99.72           N
ATOM    367  CA  CYS A  48     -41.351 -14.915   6.380  1.00106.82           C
ATOM    368  C   CYS A  48     -41.146 -13.719   7.305  1.00104.36           C
ATOM    369  O   CYS A  48     -40.407 -13.794   8.283  1.00104.56           O
ATOM    370  CB  CYS A  48     -42.835 -15.275   6.303  1.00105.79           C
ATOM    371  SG  CYS A  48     -43.707 -15.186   7.879  1.00101.00           S
TER
HETATM  459 ZN    ZN A 102     -44.322 -17.446   8.351  1.00 98.67          ZN
END
"""

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
