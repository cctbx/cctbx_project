from __future__ import absolute_import, division, print_function
import time
from mmtbx.hydrogens.tst_add_hydrogen_1 import compare_models

def run():
  test_001()
  test_009()

# ------------------------------------------------------------------------------

def test_001():
  '''
    CD1 is missing --> H, HD1, HE1 can't be parameterized --> should not be added.
  '''
  compare_models(pdb_str      = pdb_str_001,
                 not_contains = ' HE1')

# ------------------------------------------------------------------------------

def test_009():
  '''
    7JX4: ACE linked to GLY1 Goal: H1,2,3 on ACE and no H to HYP.
  '''
  compare_models(pdb_str = pdb_str_009)


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

pdb_str_009 = '''
REMARK 7JX4: ACE linked at GLY 1 and HYP (PRO-like)
CRYST1   72.365   24.756   25.357  90.00  98.72  90.00 C 1 2 1
HETATM    1  C   ACE A   0      38.724   1.027  10.647  1.00 35.12           C
HETATM    2  O   ACE A   0      39.785   1.620  10.713  1.00 31.48           O
HETATM    3  CH3 ACE A   0      37.707   1.320   9.557  1.00 37.70           C
HETATM    4  H1  ACE A   0      37.394   2.362   9.563  1.00 37.70           H
HETATM    5  H2  ACE A   0      36.807   0.718   9.668  1.00 37.70           H
HETATM    6  H3  ACE A   0      38.104   1.115   8.565  1.00 37.70           H
ATOM      7  N   GLY A   1      38.359   0.034  11.445  1.00 26.83           N
ATOM      8  CA  GLY A   1      39.087  -0.442  12.689  1.00 24.95           C
ATOM      9  C   GLY A   1      40.351   0.366  12.956  1.00 20.49           C
ATOM     10  O   GLY A   1      40.732   1.220  12.165  1.00 21.71           O
ATOM     11  H   GLY A   1      37.650  -0.432  11.306  1.00 26.83           H
ATOM     12  HA2 GLY A   1      38.501  -0.358  13.457  1.00 24.95           H
ATOM     13  HA3 GLY A   1      39.335  -1.374  12.581  1.00 24.95           H
ATOM     14  N   PRO A   2      41.104   0.126  14.061  1.00 20.40           N
ATOM     15  CA  PRO A   2      42.349   0.839  14.257  1.00 19.14           C
ATOM     16  C   PRO A   2      43.339   0.479  13.161  1.00 15.63           C
ATOM     17  O   PRO A   2      43.259  -0.533  12.439  1.00 14.90           O
ATOM     18  CB  PRO A   2      42.883   0.369  15.625  1.00 23.56           C
ATOM     19  CG  PRO A   2      42.199  -0.976  15.820  1.00 27.80           C
ATOM     20  CD  PRO A   2      40.839  -0.836  15.147  1.00 24.02           C
ATOM     21  HA  PRO A   2      42.194   1.795  14.310  1.00 19.14           H
ATOM     22  HB2 PRO A   2      43.848   0.274  15.596  1.00 23.56           H
ATOM     23  HB3 PRO A   2      42.632   0.998  16.320  1.00 23.56           H
ATOM     24  HG2 PRO A   2      42.722  -1.676  15.399  1.00 27.80           H
ATOM     25  HG3 PRO A   2      42.099  -1.160  16.767  1.00 27.80           H
ATOM     26  HD2 PRO A   2      40.180  -0.486  15.767  1.00 24.02           H
ATOM     27  HD3 PRO A   2      40.540  -1.687  14.791  1.00 24.02           H
HETATM   28  N   HYP A   3      44.381   1.316  13.052  1.00 14.81           N
HETATM   29  CA  HYP A   3      45.534   1.014  12.227  1.00 13.64           C
HETATM   30  C   HYP A   3      46.092  -0.340  12.609  1.00 11.90           C
HETATM   31  O   HYP A   3      46.074  -0.816  13.778  1.00 13.30           O
HETATM   32  CB  HYP A   3      46.514   2.141  12.599  1.00 15.74           C
HETATM   33  CG  HYP A   3      45.629   3.294  13.029  1.00 18.10           C
HETATM   34  CD  HYP A   3      44.531   2.590  13.781  1.00 18.00           C
HETATM   35  OD1 HYP A   3      45.026   3.949  11.913  1.00 18.20           O
HETATM   37  HA  HYP A   3      45.349   1.056  11.276  1.00 13.64           H
HETATM   38  HB2 HYP A   3      47.111   1.889  13.321  1.00 15.74           H
HETATM   39  HB3 HYP A   3      47.072   2.408  11.852  1.00 15.74           H
HETATM   40  HG  HYP A   3      46.135   3.928  13.562  1.00 18.10           H
HETATM   41  HD1 HYP A   3      44.639   3.368  11.428  1.00 18.20           H
HETATM   42 HD22 HYP A   3      44.765   2.430  14.709  1.00 18.00           H
HETATM   43 HD23 HYP A   3      43.700   3.091  13.773  1.00 18.00           H
ATOM     44  N   GLY A   4      46.667  -0.918  11.580  1.00 10.55           N
ATOM     45  CA  GLY A   4      47.460  -2.080  11.829  1.00 10.24           C
ATOM     46  C   GLY A   4      48.729  -1.736  12.583  1.00  9.86           C
ATOM     47  O   GLY A   4      49.036  -0.544  12.829  1.00 11.39           O
ATOM     48  OXT GLY A   4      49.480  -2.633  12.966  1.00  9.86           O
ATOM     49  H   GLY A   4      46.611  -0.661  10.761  1.00 10.55           H
ATOM     50  HA2 GLY A   4      46.949  -2.714  12.356  1.00 10.24           H
ATOM     51  HA3 GLY A   4      47.703  -2.493  10.986  1.00 10.24           H
'''


if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
