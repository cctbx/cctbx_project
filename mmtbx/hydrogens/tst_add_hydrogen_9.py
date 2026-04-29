from __future__ import absolute_import, division, print_function
import time
from mmtbx.hydrogens.tst_add_hydrogen_1 import compare_models

def run():
  test_001()

# ------------------------------------------------------------------------------

def test_001():
  '''
    CD1 is missing --> H, HD1, HE1 can't be parameterized --> should not be added.
  '''
  compare_models(pdb_str      = pdb_str_001,
                 not_contains = ' HE1')

# ------------------------------------------------------------------------------

pdb_str_001 = """
REMARK CD1 is missing --> HD1, HE1 cannot not be placed; H is not placed
CRYST1   17.955   13.272   13.095  90.00  90.00  90.00 P 1
ATOM      1  N   TYR A 139      10.241   7.920   5.000  1.00 10.00           N
ATOM      2  CA  TYR A 139      10.853   7.555   6.271  1.00 10.00           C
ATOM      3  C   TYR A 139      12.362   7.771   6.227  1.00 10.00           C
ATOM      4  O   TYR A 139      12.955   8.272   7.181  1.00 10.00           O
ATOM      5  CB  TYR A 139      10.540   6.098   6.617  1.00 10.00           C
ATOM      6  CG  TYR A 139       9.063   5.805   6.749  1.00 10.00           C
ATOM      7  CD2 TYR A 139       8.414   5.943   7.969  1.00 10.00           C
ATOM      8  CE1 TYR A 139       6.966   5.122   5.770  1.00 10.00           C
ATOM      9  CE2 TYR A 139       7.064   5.676   8.095  1.00 10.00           C
ATOM     10  CZ  TYR A 139       6.345   5.266   6.993  1.00 10.00           C
ATOM     11  OH  TYR A 139       5.000   5.000   7.113  1.00 10.00           O
ATOM     12  H   TYR A 139       9.382   7.879   5.001  1.00 10.00           H
ATOM     13  HA  TYR A 139      10.487   8.115   6.973  1.00 10.00           H
ATOM     14  HB2 TYR A 139      10.961   5.881   7.464  1.00 10.00           H
ATOM     15  HB3 TYR A 139      10.893   5.529   5.916  1.00 10.00           H
ATOM     16  HD2 TYR A 139       8.896   6.220   8.714  1.00 10.00           H
ATOM     17  HE2 TYR A 139       6.643   5.772   8.919  1.00 10.00           H
ATOM     18  HH  TYR A 139       4.752   5.127   7.905  1.00 10.00           H
TER
"""

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
