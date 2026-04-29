from __future__ import absolute_import, division, print_function
import time

from mmtbx.hydrogens.tst_add_hydrogen_1 import compare_models

# ------------------------------------------------------------------------------

def run():
  test_000()
  test_001()
  test_002()
  test_003()

# ------------------------------------------------------------------------------

def test_000():
  '''
  D-Leucine with an NH2.
  '''
  compare_models(pdb_str = pdb_str_000)

# ------------------------------------------------------------------------------

def test_001():
  '''
  NAG carbohydrate
  '''
  compare_models(pdb_str = pdb_str_001)

# ------------------------------------------------------------------------------

def test_002():
  '''
  ADM linked to GLY. Make sure HA2 is tetrahedral.
  '''
  compare_models(pdb_str = pdb_str_002)

# ------------------------------------------------------------------------------

def test_003():
  '''
  GLY1 linked to the sidechain of ASP9. Make sure only one H of the propeller
  remains and that it has correct geometry (planar).
  '''
  compare_models(pdb_str = pdb_str_003)

# ------------------------------------------------------------------------------

pdb_str_000 = """
CRYST1   56.445   72.085  123.593  90.00  90.00  90.00 P 21 21 21
SCALE1      0.017716  0.000000  0.000000        0.00000
SCALE2      0.000000  0.013873  0.000000        0.00000
SCALE3      0.000000  0.000000  0.008091        0.00000
HETATM    1  N   DLE E  12     -46.710 -11.701  15.468  1.00 34.07           N
HETATM    2  CA  DLE E  12     -47.200 -11.850  14.105  1.00 39.64           C
HETATM    3  C   DLE E  12     -48.679 -11.482  14.009  1.00 41.13           C
HETATM    4  O   DLE E  12     -49.335 -11.798  13.021  1.00 42.80           O
HETATM    5  CB  DLE E  12     -46.390 -10.976  13.144  1.00 43.98           C
HETATM    6  CG  DLE E  12     -44.950 -11.415  12.889  1.00 46.65           C
HETATM    7  CD1 DLE E  12     -44.919 -12.765  12.197  1.00 46.88           C
HETATM    8  CD2 DLE E  12     -44.198 -10.363  12.085  1.00 48.25           C
HETATM    9  H   DLE E  12     -46.283 -12.387  15.763  1.00 34.07           H
HETATM   10  HA  DLE E  12     -47.098 -12.777  13.837  1.00 39.64           H
HETATM   11  HB2 DLE E  12     -46.356 -10.078  13.509  1.00 43.98           H
HETATM   12  HB3 DLE E  12     -46.844 -10.969  12.287  1.00 43.98           H
HETATM   13  HG  DLE E  12     -44.493 -11.511  13.739  1.00 46.65           H
HETATM   14 HD11 DLE E  12     -43.996 -13.021  12.046  1.00 46.88           H
HETATM   15 HD12 DLE E  12     -45.386 -12.696  11.349  1.00 46.88           H
HETATM   16 HD13 DLE E  12     -45.356 -13.420  12.763  1.00 46.88           H
HETATM   17 HD21 DLE E  12     -44.646 -10.236  11.234  1.00 48.25           H
HETATM   18 HD22 DLE E  12     -43.289 -10.668  11.938  1.00 48.25           H
HETATM   19 HD23 DLE E  12     -44.192  -9.530  12.582  1.00 48.25           H
HETATM   20  N   NH2 E 202     -49.193 -10.809  15.033  1.00 40.90           N
HETATM   21  HN1 NH2 E 202     -48.680 -10.596  15.728  1.00 40.90           H
HETATM   22  HN2 NH2 E 202     -50.051 -10.572  15.024  1.00 40.90           H
END
"""

pdb_str_001 = '''
CRYST1   16.163   16.054   17.729  90.00  90.00  90.00 P 1
SCALE1      0.061870  0.000000  0.000000        0.00000
SCALE2      0.000000  0.062290  0.000000        0.00000
SCALE3      0.000000  0.000000  0.056405        0.00000
HETATM    1  C1  NAG A   1       7.207   7.892   8.696  1.00 20.00      A    C
HETATM    2  C2  NAG A   1       8.726   7.903   8.675  1.00 20.00      A    C
HETATM    3  C3  NAG A   1       9.299   7.873  10.037  1.00 20.00      A    C
HETATM    4  C4  NAG A   1       8.735   8.963  10.912  1.00 20.00      A    C
HETATM    5  C5  NAG A   1       7.211   8.952  10.928  1.00 20.00      A    C
HETATM    6  C6  NAG A   1       6.722  10.156  11.692  1.00 20.00      A    C
HETATM    7  C7  NAG A   1       9.920   6.866   6.642  1.00 20.00      A    C
HETATM    8  C8  NAG A   1      10.391   5.633   5.851  1.00 20.00      A    C
HETATM    9  N2  NAG A   1       9.210   6.694   7.918  1.00 20.00      A    N
HETATM   10  O1  NAG A   1       6.748   8.064   7.430  1.00 20.00      A    O
HETATM   11  O3  NAG A   1      10.730   8.045   9.943  1.00 20.00      A    O
HETATM   12  O4  NAG A   1       9.211   8.778  12.243  1.00 20.00      A    O
HETATM   13  O5  NAG A   1       6.644   8.983   9.574  1.00 20.00      A    O
HETATM   14  O6  NAG A   1       5.328  10.234  11.597  1.00 20.00      A    O
HETATM   15  O7  NAG A   1      10.109   7.947   6.214  1.00 20.00      A    O
HETATM   16  H1  NAG A   1       6.928   7.042   9.072  1.00 20.00      A    H
HETATM   17  H2  NAG A   1       9.004   8.727   8.246  1.00 20.00      A    H
HETATM   18  H3  NAG A   1       9.042   7.013  10.404  1.00 20.00      A    H
HETATM   19  H4  NAG A   1       9.054   9.784  10.505  1.00 20.00      A    H
HETATM   20  H5  NAG A   1       6.893   8.131  11.334  1.00 20.00      A    H
HETATM   21  H61 NAG A   1       7.154  10.946  11.330  1.00 20.00      A    H
HETATM   22  H62 NAG A   1       7.018  10.077  12.612  1.00 20.00      A    H
HETATM   23  H81 NAG A   1      10.095   4.828   6.305  1.00 20.00      A    H
HETATM   24  H82 NAG A   1      11.359   5.638   5.795  1.00 20.00      A    H
HETATM   25  H83 NAG A   1      10.013   5.661   4.958  1.00 20.00      A    H
HETATM   26  HN2 NAG A   1       9.068   5.908   8.238  1.00 20.00      A    H
HETATM   27  HO1 NAG A   1       6.806   8.885   7.217  1.00 20.00      A    H
HETATM   28  HO3 NAG A   1      10.879   8.808   9.599  1.00 20.00      A    H
HETATM   29  HO4 NAG A   1       9.998   8.457  12.213  1.00 20.00      A    H
HETATM   30  HO6 NAG A   1       5.116  10.219  10.774  1.00 20.00      A    H
END
'''

pdb_str_002 = '''
REMARK Check that H near link has correct geometry
CRYST1   30.746   56.560   99.515  90.00  90.00  90.00 P 21 21 21
SCALE1      0.032525  0.000000  0.000000        0.00000
SCALE2      0.000000  0.017680  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010049        0.00000
ATOM      1  N   GLY B  10       2.957 -10.749  -1.364  1.00 22.36           N
ATOM      2  CA  GLY B  10       2.830 -10.787  -2.810  1.00 23.42           C
ATOM      3  C   GLY B  10       3.228 -12.114  -3.425  1.00 21.76           C
ATOM      4  O   GLY B  10       2.672 -12.516  -4.453  1.00 20.78           O
ATOM      5  H   GLY B  10       2.574 -11.405  -0.962  1.00 22.36           H
ATOM      6  HA2 GLY B  10       1.884 -10.688  -3.000  1.00 23.42           H
TER
HETATM    8  C1  ADM B 101       3.506  -9.585  -3.485  1.00 29.99           C
HETATM    9  C10 ADM B 101       3.272  -6.994  -4.853  1.00 33.86           C
HETATM   10  C2  ADM B 101       2.997  -9.452  -4.919  1.00 37.40           C
HETATM   11  C3  ADM B 101       3.647  -8.260  -5.608  1.00 30.02           C
HETATM   12  C4  ADM B 101       5.159  -8.444  -5.610  1.00 40.96           C
HETATM   13  C5  ADM B 101       5.668  -8.541  -4.180  1.00 34.29           C
HETATM   14  C6  ADM B 101       5.301  -7.254  -3.457  1.00 29.01           C
HETATM   15  C7  ADM B 101       3.789  -7.091  -3.424  1.00 35.16           C
HETATM   16  C8  ADM B 101       3.171  -8.301  -2.730  1.00 32.12           C
HETATM   17  C9  ADM B 101       5.021  -9.747  -3.503  1.00 25.37           C
HETATM   19  H3  ADM B 101       3.336  -8.190  -6.524  1.00 30.02           H
HETATM   20  H21 ADM B 101       2.033  -9.347  -4.917  1.00 37.40           H
HETATM   21  H22 ADM B 101       3.189 -10.264  -5.413  1.00 37.40           H
HETATM   22  H41 ADM B 101       5.585  -7.701  -6.065  1.00 40.96           H
HETATM   23  H42 ADM B 101       5.397  -9.245  -6.103  1.00 40.96           H
HETATM   24  H5  ADM B 101       6.631  -8.660  -4.160  1.00 34.29           H
HETATM   25  H61 ADM B 101       5.655  -7.268  -2.554  1.00 29.01           H
HETATM   26  H62 ADM B 101       5.709  -6.495  -3.902  1.00 29.01           H
HETATM   27  H7  ADM B 101       3.551  -6.285  -2.940  1.00 35.16           H
HETATM   28  H81 ADM B 101       3.497  -8.358  -1.818  1.00 32.12           H
HETATM   29  H82 ADM B 101       2.209  -8.191  -2.674  1.00 32.12           H
HETATM   30  H91 ADM B 101       5.357  -9.839  -2.598  1.00 25.37           H
HETATM   31  H92 ADM B 101       5.263 -10.560  -3.972  1.00 25.37           H
HETATM   32 H101 ADM B 101       3.648  -6.216  -5.294  1.00 33.86           H
HETATM   33 H102 ADM B 101       2.309  -6.875  -4.854  1.00 33.86           H
END
'''

pdb_str_003 = '''
REMARK based on PDB model 3njw
CRYST1   19.465   21.432   29.523  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   GLY A   1       6.011  23.726   5.538  1.00  4.36           N
ATOM      2  CA  GLY A   1       7.279  24.418   5.504  1.00  4.79           C
ATOM      3  C   GLY A   1       8.370  23.751   6.291  1.00  4.66           C
ATOM      4  O   GLY A   1       9.449  24.344   6.484  1.00  5.69           O
ATOM      5  H1  GLY A   1       5.585  23.650   6.316  1.00  4.36           H
ATOM      8  HA2 GLY A   1       7.158  25.311   5.862  1.00  4.79           H
ATOM      9  HA3 GLY A   1       7.575  24.484   4.583  1.00  4.79           H
ATOM     10  N   ASP A   9       3.861  20.962   2.856  1.00  4.41           N
ATOM     11  CA  ASP A   9       3.324  22.217   3.387  1.00  4.26           C
ATOM     12  C   ASP A   9       3.240  23.233   2.275  1.00  4.25           C
ATOM     13  O   ASP A   9       3.655  23.005   1.142  1.00  4.88           O
ATOM     14  CB  ASP A   9       4.061  22.659   4.637  1.00  3.84           C
ATOM     15  CG  ASP A   9       5.454  23.195   4.444  1.00  3.99           C
ATOM     16  OD1 ASP A   9       6.011  23.127   3.350  1.00  4.95           O
ATOM     17  H   ASP A   9       3.281  20.481   2.442  1.00  4.41           H
ATOM     18  HA  ASP A   9       2.416  22.112   3.712  1.00  4.26           H
ATOM     19  HB2 ASP A   9       3.543  23.362   5.059  1.00  3.84           H
ATOM     20  HB3 ASP A   9       4.131  21.896   5.231  1.00  3.84           H
END
'''

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
