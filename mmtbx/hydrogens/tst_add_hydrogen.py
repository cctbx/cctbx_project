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

# ------------------------------------------------------------------------------
  # make sure that H of incomplete peptide unit (N-terminal) is not placed
  # test will need to be adapted once H3 can be placed at N-terminal
# ------------------------------------------------------------------------------
def exercise1():
  pdb_inp = iotbx.pdb.input(lines=pdb_str_1.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())

  model_with_h = mmtbx.hydrogens.add(model = model)

  number_h = model_with_h.get_hd_selection().count(True)
  # total number of H atoms in Tyr: 9
  # but peptide H atom cannot be parameterized
  assert(number_h == 8)

# ------------------------------------------------------------------------------
  # if one heavy atom is missing, it can impact parameterization and thus
  # placement of SEVERAL H atoms.
# ------------------------------------------------------------------------------
def exercise2():
  pdb_inp = iotbx.pdb.input(lines=pdb_str_2.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp,
                              log=null_out())

  model_with_h = mmtbx.hydrogens.add(model = model)

  number_h = model_with_h.get_hd_selection().count(True)
  # The following atoms cannot be parameterized:
  # H, HD1, HE1
  # number_h = 9 - 3 = 6
  assert(number_h == 6)

  geometry_restraints = model_with_h.get_restraints_manager().geometry
  e_sites = geometry_restraints.energies_sites(
          sites_cart        = model_with_h.get_sites_cart(),
          compute_gradients = False)
  t = e_sites.target
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

#  print(e_sites.bond_residual_sum)
#  print(e_sites.angle_residual_sum)
#  print(e_sites.residual_sum)

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

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
