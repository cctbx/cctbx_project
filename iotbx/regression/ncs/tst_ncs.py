from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal, show_diff
from scitbx import matrix
from iotbx import pdb
import mmtbx.model
import iotbx.ncs
import iotbx.ncs as ncs

pdb_str_1="""\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
ATOM      1  N   THR A   1       9.670  10.289  11.135  1.00 20.00           N
ATOM      2  CA  THR A   2       9.559   8.931  10.615  1.00 20.00           C
ATOM      3  C   THR A   3       9.634   7.903  11.739  1.00 20.00           C
ATOM      4  O   THR B   4      10.449   8.027  12.653  1.00 20.00           O
ATOM      5  CB  THR B   5      10.660   8.630   9.582  1.00 20.00           C
ATOM      6  OG1 THR A   6      10.560   9.552   8.490  1.00 20.00           O
ATOM      7  CG2 THR A   7      10.523   7.209   9.055  1.00 20.00           C
TER
"""

pdb_str_3="""
REMARK   0 Test molecule with BIOMOLECULE: 1
REMARK   0
REMARK   0 The test will generate the biomolecule (the multimer assembly)
REMARK   0 from the transformation matrices writen below
REMARK   0 and then compare the results to the calculated expected one
REMARK 350 CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN.
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.000000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.000000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT1   2  1.000000  0.000000  0.000000        0.000000
REMARK 350   BIOMT2   2  0.000000  0.000000 -1.000000        0.000000
REMARK 350   BIOMT3   2  0.000000  1.000000  0.000000        0.000000
REMARK 350   BIOMT1   3  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT2   3  0.000000  1.000000  0.000000        0.000000
REMARK 350   BIOMT3   3 -1.000000  0.000000  0.000000        0.000000
REMARK 350   BIOMT1   4  0.000000 -1.000000  0.000000        0.000000
REMARK 350   BIOMT2   4  1.000000  0.000000  0.000000        0.000000
REMARK 350   BIOMT3   4  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT1   5  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT2   5  0.000000  1.000000  0.000000        0.000000
REMARK 350   BIOMT3   5 -1.000000  0.000000  0.000000        0.000000
REMARK 350   BIOMT1   6  0.000000 -1.000000  0.000000        0.000000
REMARK 350   BIOMT2   6  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT3   6 -1.000000  0.000000  0.000000        0.000000
REMARK 350   BIOMT1   7  0.500000 -0.866025  0.000000        0.000000
REMARK 350   BIOMT2   7  0.866025  0.500000  0.000000        0.000000
REMARK 350   BIOMT3   7  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT1   8 -0.500000 -0.866025  0.000000        0.000000
REMARK 350   BIOMT2   8  0.866025 -0.500000  0.000000        0.000000
REMARK 350   BIOMT3   8  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT1   9  1.000000  0.000000  0.000000        0.000000
REMARK 350   BIOMT2   9  0.000000  1.000000  0.000000        0.500000
REMARK 350   BIOMT3   9  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT1  10 -0.500000 -0.866025  0.000000        0.000000
REMARK 350   BIOMT2  10  0.866025 -0.500000  0.000000        0.000000
REMARK 350   BIOMT3  10  0.000000  0.000000  1.000000       -1.000000
MTRIX1   1  1.000000  0.000000  0.000000        0.00000
MTRIX2   1  0.000000  1.000000  0.000000        0.00000
MTRIX3   1  0.000000  0.000000  1.000000        0.00000
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  1.000000  0.000000  0.000000        0.00000
MTRIX2   2  0.000000  0.000000 -1.000000        0.00000
MTRIX3   2  0.000000  1.000000  0.000000        0.00000
MTRIX1   3  0.500000 -0.866025  0.000000        0.00000
MTRIX2   3  0.866025  0.500000  0.000000        0.00000
MTRIX3   3  0.000000  0.000000  1.000000        0.00000
MTRIX1   4 -0.500000 -0.866025  0.000000        0.00000
MTRIX2   4  0.866025 -0.500000  0.000000        0.00000
MTRIX3   4  0.000000  0.000000  1.000000        0.00000
MTRIX1   5  1.000000  0.000000  0.000000        0.00000
MTRIX2   5  0.000000  1.000000  0.000000        0.50000
MTRIX3   5  0.000000  0.000000  1.000000        0.00000
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
ATOM      1   N  ILE A  40       1.000   1.000   1.000  1.00162.33           C
ATOM      2  CA  LEU A  40      94.618  -5.253  91.582  1.00 87.10           C
ATOM      3   C  ARG B  40      62.395  51.344  80.786  1.00107.25           C
HETATM    4  C1  EDO A  40      39.954  51.526  72.372  0.33 60.93           C
"""

pdb_str_4 = """\
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2  0.309017 -0.951057  0.000000        0.00000
REMARK 350   BIOMT2   2  0.951057  0.309017 -0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        7.00000
REMARK 350   BIOMT1   3 -0.809017 -0.587785  0.000000        0.00000
REMARK 350   BIOMT2   3  0.587785 -0.809017 -0.000000        0.00000
REMARK 350   BIOMT3   3  0.000000  0.000000  1.000000        0.00000
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   2      64.807-112.186 260.746  1.00160.99           N
ATOM      2  CA  ALA A   2      64.727-111.450 262.002  1.00159.36           C
ATOM      3  C   ALA A   2      63.960-110.148 261.805  1.00154.38           C
ATOM      4  O   ALA A   2      62.935-109.914 262.452  1.00149.47           O
ATOM      5  CB  ALA A   2      66.123-111.175 262.542  1.00156.98           C
ATOM      6  N   SER A   3      64.474-109.323 260.896  1.00135.75           N
ATOM      7  CA  SER A   3      63.887-108.040 260.510  1.00131.97           C
ATOM      8  C   SER A   3      64.863-107.340 259.575  1.00140.51           C
ATOM      9  O   SER A   3      65.864-107.925 259.165  1.00148.46           O
ATOM     10  CB  SER A   3      63.641-107.147 261.726  1.00126.01           C
ATOM     11  OG  SER A   3      64.002-105.804 261.453  1.00119.04           O
END
"""

pdb_str_5 = """\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.309017 -0.951057  0.000000        0.00000
MTRIX2   2  0.951057  0.309017 -0.000000        0.00000
MTRIX3   2  0.000000  0.000000  1.000000        7.00000
MTRIX1   3 -0.809017 -0.587785  0.000000        0.00000
MTRIX2   3  0.587785 -0.809017 -0.000000        0.00000
MTRIX3   3  0.000000  0.000000  1.000000        0.00000
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
ATOM    757  N   ASP A 247      16.068 -20.882 -28.984  1.00 35.93           N
ATOM    758  CA  ASP A 247      15.914 -22.265 -28.600  1.00 47.90           C
ATOM    759  C   ASP A 247      17.130 -23.042 -29.116  1.00 42.32           C
ATOM    760  O   ASP A 247      17.461 -22.986 -30.301  1.00 47.25           O
ATOM    761  CB  ASP A 247      14.621 -22.814 -29.198  1.00 47.22           C
ATOM    762  CG  ASP A 247      14.068 -23.974 -28.412  1.00 61.15           C
ATOM    763  OD1 ASP A 247      14.359 -24.061 -27.196  1.00 63.66           O
ATOM    764  OD2 ASP A 247      13.341 -24.798 -29.012  1.00 77.01           O
ATOM    765  N   VAL A 248      17.808 -23.746 -28.218  1.00 44.08           N
ATOM    766  CA  VAL A 248      19.008 -24.503 -28.584  1.00 46.18           C
ATOM    767  C   VAL A 248      18.668 -25.988 -28.583  1.00 53.97           C
ATOM    768  O   VAL A 248      18.049 -26.478 -27.638  1.00 51.48           O
ATOM    769  CB  VAL A 248      20.185 -24.226 -27.608  1.00 47.55           C
ATOM    770  CG1 VAL A 248      21.414 -25.015 -28.012  1.00 41.43           C
ATOM    771  CG2 VAL A 248      20.513 -22.743 -27.567  1.00 41.64           C
ATOM    772  N   VAL A 249      19.057 -26.697 -29.641  1.00 54.29           N
ATOM    773  CA  VAL A 249      18.662 -28.097 -29.810  1.00 60.17           C
ATOM    774  C   VAL A 249      19.859 -29.041 -29.982  1.00 57.98           C
ATOM    775  O   VAL A 249      20.731 -28.827 -30.828  1.00 58.31           O
ATOM    776  CB  VAL A 249      17.671 -28.280 -30.997  1.00 60.85           C
ATOM    777  CG1 VAL A 249      16.500 -27.300 -30.884  1.00 48.00           C
ATOM    778  CG2 VAL A 249      18.386 -28.110 -32.337  1.00 59.99           C
TER
ATOM    780  N   LYS B 151       4.045  -6.858 -32.823  1.00 45.22           N
ATOM    781  CA  LYS B 151       4.686  -6.715 -34.123  1.00 50.40           C
ATOM    782  C   LYS B 151       5.707  -5.554 -34.172  1.00 47.13           C
ATOM    783  O   LYS B 151       6.820  -5.764 -34.625  1.00 52.91           O
ATOM    784  CB  LYS B 151       3.657  -6.646 -35.268  1.00 40.73           C
ATOM    785  CG  LYS B 151       4.264  -6.627 -36.661  1.00 55.98           C
ATOM    786  CD  LYS B 151       3.272  -7.051 -37.745  1.00 72.14           C
ATOM    787  CE  LYS B 151       2.529  -8.338 -37.375  1.00 75.11           C
ATOM    788  NZ  LYS B 151       3.451  -9.400 -36.884  1.00 75.46           N
ATOM    789  N   ARG B 152       5.369  -4.349 -33.709  1.00 42.01           N
ATOM    790  CA  ARG B 152       6.399  -3.290 -33.702  1.00 40.51           C
ATOM    791  C   ARG B 152       6.155  -2.002 -32.909  1.00 34.21           C
ATOM    792  O   ARG B 152       5.015  -1.605 -32.636  1.00 33.77           O
ATOM    793  CB  ARG B 152       6.845  -2.937 -35.130  1.00 40.62           C
ATOM    794  CG  ARG B 152       5.842  -2.126 -35.925  1.00 45.94           C
ATOM    795  CD  ARG B 152       6.341  -1.926 -37.341  1.00 42.75           C
ATOM    796  NE  ARG B 152       7.478  -1.006 -37.404  1.00 45.27           N
ATOM    797  CZ  ARG B 152       8.177  -0.763 -38.509  1.00 49.68           C
ATOM    798  NH1 ARG B 152       7.860  -1.382 -39.644  1.00 47.81           N
ATOM    799  NH2 ARG B 152       9.192   0.096 -38.482  1.00 48.06           N
END
"""

pdb_str_8 = """\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000    1
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000    1
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000    1
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000    1
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000    1
MTRIX3   3  0.565855  0.754120  0.333333        0.00000    1
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
ATOM      1  N   THR A   1       9.670  10.289  11.135  1.00 20.00           N
ATOM      2  CA  THR A   2       9.559   8.931  10.615  1.00 20.00           C
ATOM      3  C   THR A   3       9.634   7.903  11.739  1.00 20.00           C
ATOM      4  O   THR B   4      10.449   8.027  12.653  1.00 20.00           O
ATOM      5  CB  THR B   5      10.660   8.630   9.582  1.00 20.00           C
ATOM      6  OG1 THR A   6      10.560   9.552   8.490  1.00 20.00           O
ATOM      7  CG2 THR A   7      10.523   7.209   9.055  1.00 20.00           C
END
"""

# A dimer
pdb_str_9=\
"""
CRYST1   18.411   24.923   19.598  90.00  90.00  90.00 P 1
ATOM      1  N   LEU A 219      12.476   8.570   6.342  1.00 30.00           N
ATOM      2  CA  LEU A 219      11.029   8.720   6.198  1.00 30.00           C
ATOM      3  C   LEU A 219      10.371   9.171   7.496  1.00 30.00           C
ATOM      4  O   LEU A 219       9.422   9.957   7.469  1.00 30.00           O
ATOM      5  CB  LEU A 219      10.404   7.413   5.700  1.00 30.00           C
ATOM      6  CG  LEU A 219       8.894   7.444   5.410  1.00 30.00           C
ATOM      7  CD1 LEU A 219       8.567   6.777   4.078  1.00 30.00           C
ATOM      8  CD2 LEU A 219       8.092   6.808   6.540  1.00 30.00           C
ATOM      9  N   ARG A 220      10.878   8.672   8.622  1.00 30.00           N
ATOM     10  CA  ARG A 220      10.381   9.053   9.938  1.00 30.00           C
ATOM     11  C   ARG A 220      10.603  10.531  10.185  1.00 30.00           C
ATOM     12  O   ARG A 220       9.715  11.202  10.728  1.00 30.00           O
ATOM     13  CB  ARG A 220      11.065   8.236  11.041  1.00 30.00           C
ATOM     14  CG  ARG A 220      10.661   6.771  11.083  1.00 30.00           C
ATOM     15  CD  ARG A 220       9.321   6.571  11.776  1.00 30.00           C
ATOM     16  NE  ARG A 220       8.945   5.157  11.847  1.00 30.00           N
ATOM     17  CZ  ARG A 220       9.446   4.259  12.701  1.00 30.00           C
ATOM     18  NH1 ARG A 220      10.377   4.593  13.598  1.00 30.00           N
ATOM     19  NH2 ARG A 220       9.012   3.000  12.656  1.00 30.00           N
ATOM     20  N   LEU A 221      11.792  11.006   9.803  1.00 30.00           N
ATOM     21  CA  LEU A 221      12.179  12.402   9.990  1.00 30.00           C
ATOM     22  C   LEU A 221      11.192  13.323   9.293  1.00 30.00           C
ATOM     23  O   LEU A 221      10.739  14.325   9.874  1.00 30.00           O
ATOM     24  CB  LEU A 221      13.578  12.678   9.432  1.00 30.00           C
ATOM     25  CG  LEU A 221      14.169  14.072   9.673  1.00 30.00           C
ATOM     26  CD1 LEU A 221      14.503  14.277  11.143  1.00 30.00           C
ATOM     27  CD2 LEU A 221      15.411  14.272   8.823  1.00 30.00           C
ATOM     28  N   GLN A 222      10.879  12.956   8.053  1.00 30.00           N
ATOM     29  CA  GLN A 222       9.943  13.723   7.218  1.00 30.00           C
ATOM     30  C   GLN A 222       8.591  13.825   7.903  1.00 30.00           C
ATOM     31  O   GLN A 222       7.994  14.910   7.968  1.00 30.00           O
ATOM     32  CB  GLN A 222       9.767  13.081   5.840  1.00 30.00           C
ATOM     33  CG  GLN A 222       8.924  13.893   4.867  1.00 30.00           C
ATOM     34  CD  GLN A 222       8.715  13.178   3.545  1.00 30.00           C
ATOM     35  OE1 GLN A 222       9.626  12.537   3.018  1.00 30.00           O
ATOM     36  NE2 GLN A 222       7.507  13.285   3.000  1.00 30.00           N
ATOM     37  N   GLN A 223       8.143  12.676   8.406  1.00 30.00           N
ATOM     38  CA  GLN A 223       6.857  12.566   9.106  1.00 30.00           C
ATOM     39  C   GLN A 223       6.827  13.513  10.297  1.00 30.00           C
ATOM     40  O   GLN A 223       5.844  14.249  10.510  1.00 30.00           O
ATOM     41  CB  GLN A 223       6.617  11.147   9.642  1.00 30.00           C
ATOM     42  CG  GLN A 223       6.137  10.118   8.636  1.00 30.00           C
ATOM     43  CD  GLN A 223       5.812   8.791   9.300  1.00 30.00           C
ATOM     44  OE1 GLN A 223       6.419   7.765   8.994  1.00 30.00           O
ATOM     45  NE2 GLN A 223       4.859   8.809  10.228  1.00 30.00           N
ATOM     46  N   ALA A 224       7.918  13.468  11.053  1.00 30.00           N
ATOM     47  CA  ALA A 224       8.090  14.293  12.253  1.00 30.00           C
ATOM     48  C   ALA A 224       7.967  15.770  11.893  1.00 30.00           C
ATOM     49  O   ALA A 224       7.262  16.539  12.568  1.00 30.00           O
ATOM     50  CB  ALA A 224       9.437  14.024  12.914  1.00 30.00           C
ATOM     51  N   TYR A 225       8.666  16.128  10.822  1.00 30.00           N
ATOM     52  CA  TYR A 225       8.683  17.504  10.316  1.00 30.00           C
ATOM     53  C   TYR A 225       7.279  17.961   9.973  1.00 30.00           C
ATOM     54  O   TYR A 225       6.865  19.071  10.346  1.00 30.00           O
ATOM     55  CB  TYR A 225       9.677  17.738   9.178  1.00 30.00           C
ATOM     56  CG  TYR A 225      10.896  18.350   9.759  1.00 30.00           C
ATOM     57  CD1 TYR A 225      10.898  19.695  10.102  1.00 30.00           C
ATOM     58  CD2 TYR A 225      12.010  17.578  10.072  1.00 30.00           C
ATOM     59  CE1 TYR A 225      11.998  20.273  10.696  1.00 30.00           C
ATOM     60  CE2 TYR A 225      13.126  18.149  10.660  1.00 30.00           C
ATOM     61  CZ  TYR A 225      13.112  19.496  10.970  1.00 30.00           C
ATOM     62  OH  TYR A 225      14.197  20.075  11.563  1.00 30.00           O
ATOM     63  N   LEU A 226       6.569  17.081   9.279  1.00 30.00           N
ATOM     64  CA  LEU A 226       5.188  17.332   8.848  1.00 30.00           C
ATOM     65  C   LEU A 226       4.312  17.613  10.059  1.00 30.00           C
ATOM     66  O   LEU A 226       3.523  18.584  10.065  1.00 30.00           O
ATOM     67  CB  LEU A 226       4.587  16.155   8.061  1.00 30.00           C
ATOM     68  CG  LEU A 226       4.521  16.312   6.542  1.00 30.00           C
ATOM     69  CD1 LEU A 226       5.914  16.258   5.943  1.00 30.00           C
ATOM     70  CD2 LEU A 226       3.643  15.232   5.923  1.00 30.00           C
ATOM     71  N   ILE A 227       4.476  16.749  11.064  1.00 30.00           N
ATOM     72  CA  ILE A 227       3.690  16.855  12.301  1.00 30.00           C
ATOM     73  C   ILE A 227       3.952  18.198  12.969  1.00 30.00           C
ATOM     74  O   ILE A 227       3.000  18.870  13.405  1.00 30.00           O
ATOM     75  CB  ILE A 227       3.668  15.644  13.279  1.00 30.00           C
ATOM     76  CG1 ILE A 227       5.016  15.359  13.918  1.00 30.00           C
ATOM     77  CG2 ILE A 227       3.107  14.410  12.588  1.00 30.00           C
ATOM     78  CD1 ILE A 227       4.926  14.502  15.160  1.00 30.00           C
ATOM     79  N   MET A 228       5.227  18.566  13.017  1.00 30.00           N
ATOM     80  CA  MET A 228       5.663  19.830  13.617  1.00 30.00           C
ATOM     81  C   MET A 228       4.991  21.002  12.916  1.00 30.00           C
ATOM     82  O   MET A 228       4.472  21.923  13.572  1.00 30.00           O
ATOM     83  CB  MET A 228       7.192  19.985  13.595  1.00 30.00           C
ATOM     84  CG  MET A 228       7.960  19.118  14.591  1.00 30.00           C
ATOM     85  SD  MET A 228       7.511  19.311  16.335  1.00 30.00           S
ATOM     86  CE  MET A 228       7.660  21.082  16.598  1.00 30.00           C
TER

ATOM      1  N   LEU B 219     -12.476  -8.570   6.342  1.00 30.00           N
ATOM      2  CA  LEU B 219     -11.029  -8.720   6.198  1.00 30.00           C
ATOM      3  C   LEU B 219     -10.371  -9.171   7.496  1.00 30.00           C
ATOM      4  O   LEU B 219      -9.422  -9.957   7.469  1.00 30.00           O
ATOM      5  CB  LEU B 219     -10.404  -7.413   5.700  1.00 30.00           C
ATOM      6  CG  LEU B 219      -8.894  -7.444   5.410  1.00 30.00           C
ATOM      7  CD1 LEU B 219      -8.567  -6.777   4.078  1.00 30.00           C
ATOM      8  CD2 LEU B 219      -8.092  -6.808   6.540  1.00 30.00           C
ATOM      9  N   ARG B 220     -10.878  -8.672   8.622  1.00 30.00           N
ATOM     10  CA  ARG B 220     -10.381  -9.053   9.938  1.00 30.00           C
ATOM     11  C   ARG B 220     -10.603 -10.531  10.185  1.00 30.00           C
ATOM     12  O   ARG B 220      -9.715 -11.202  10.728  1.00 30.00           O
ATOM     13  CB  ARG B 220     -11.065  -8.236  11.041  1.00 30.00           C
ATOM     14  CG  ARG B 220     -10.661  -6.771  11.083  1.00 30.00           C
ATOM     15  CD  ARG B 220      -9.321  -6.571  11.776  1.00 30.00           C
ATOM     16  NE  ARG B 220      -8.945  -5.157  11.847  1.00 30.00           N
ATOM     17  CZ  ARG B 220      -9.446  -4.259  12.701  1.00 30.00           C
ATOM     18  NH1 ARG B 220     -10.377  -4.593  13.598  1.00 30.00           N
ATOM     19  NH2 ARG B 220      -9.012  -3.000  12.656  1.00 30.00           N
ATOM     20  N   LEU B 221     -11.792 -11.006   9.803  1.00 30.00           N
ATOM     21  CA  LEU B 221     -12.179 -12.402   9.990  1.00 30.00           C
ATOM     22  C   LEU B 221     -11.192 -13.323   9.293  1.00 30.00           C
ATOM     23  O   LEU B 221     -10.739 -14.325   9.874  1.00 30.00           O
ATOM     24  CB  LEU B 221     -13.578 -12.678   9.432  1.00 30.00           C
ATOM     25  CG  LEU B 221     -14.169 -14.072   9.673  1.00 30.00           C
ATOM     26  CD1 LEU B 221     -14.503 -14.277  11.143  1.00 30.00           C
ATOM     27  CD2 LEU B 221     -15.411 -14.272   8.823  1.00 30.00           C
ATOM     28  N   GLN B 222     -10.879 -12.956   8.053  1.00 30.00           N
ATOM     29  CA  GLN B 222      -9.943 -13.723   7.218  1.00 30.00           C
ATOM     30  C   GLN B 222      -8.591 -13.825   7.903  1.00 30.00           C
ATOM     31  O   GLN B 222      -7.994 -14.910   7.968  1.00 30.00           O
ATOM     32  CB  GLN B 222      -9.767 -13.081   5.840  1.00 30.00           C
ATOM     33  CG  GLN B 222      -8.924 -13.893   4.867  1.00 30.00           C
ATOM     34  CD  GLN B 222      -8.715 -13.178   3.545  1.00 30.00           C
ATOM     35  OE1 GLN B 222      -9.626 -12.537   3.018  1.00 30.00           O
ATOM     36  NE2 GLN B 222      -7.507 -13.285   3.000  1.00 30.00           N
ATOM     37  N   GLN B 223      -8.143 -12.676   8.406  1.00 30.00           N
ATOM     38  CA  GLN B 223      -6.857 -12.566   9.106  1.00 30.00           C
ATOM     39  C   GLN B 223      -6.827 -13.513  10.297  1.00 30.00           C
ATOM     40  O   GLN B 223      -5.844 -14.249  10.510  1.00 30.00           O
ATOM     41  CB  GLN B 223      -6.617 -11.147   9.642  1.00 30.00           C
ATOM     42  CG  GLN B 223      -6.137 -10.118   8.636  1.00 30.00           C
ATOM     43  CD  GLN B 223      -5.812  -8.791   9.300  1.00 30.00           C
ATOM     44  OE1 GLN B 223      -6.419  -7.765   8.994  1.00 30.00           O
ATOM     45  NE2 GLN B 223      -4.859  -8.809  10.228  1.00 30.00           N
ATOM     46  N   ALA B 224      -7.918 -13.468  11.053  1.00 30.00           N
ATOM     47  CA  ALA B 224      -8.090 -14.293  12.253  1.00 30.00           C
ATOM     48  C   ALA B 224      -7.967 -15.770  11.893  1.00 30.00           C
ATOM     49  O   ALA B 224      -7.262 -16.539  12.568  1.00 30.00           O
ATOM     50  CB  ALA B 224      -9.437 -14.024  12.914  1.00 30.00           C
ATOM     51  N   TYR B 225      -8.666 -16.128  10.822  1.00 30.00           N
ATOM     52  CA  TYR B 225      -8.683 -17.504  10.316  1.00 30.00           C
ATOM     53  C   TYR B 225      -7.279 -17.961   9.973  1.00 30.00           C
ATOM     54  O   TYR B 225      -6.865 -19.071  10.346  1.00 30.00           O
ATOM     55  CB  TYR B 225      -9.677 -17.738   9.178  1.00 30.00           C
ATOM     56  CG  TYR B 225     -10.896 -18.350   9.759  1.00 30.00           C
ATOM     57  CD1 TYR B 225     -10.898 -19.695  10.102  1.00 30.00           C
ATOM     58  CD2 TYR B 225     -12.010 -17.578  10.072  1.00 30.00           C
ATOM     59  CE1 TYR B 225     -11.998 -20.273  10.696  1.00 30.00           C
ATOM     60  CE2 TYR B 225     -13.126 -18.149  10.660  1.00 30.00           C
ATOM     61  CZ  TYR B 225     -13.112 -19.496  10.970  1.00 30.00           C
ATOM     62  OH  TYR B 225     -14.197 -20.075  11.563  1.00 30.00           O
ATOM     63  N   LEU B 226      -6.569 -17.081   9.279  1.00 30.00           N
ATOM     64  CA  LEU B 226      -5.188 -17.332   8.848  1.00 30.00           C
ATOM     65  C   LEU B 226      -4.312 -17.613  10.059  1.00 30.00           C
ATOM     66  O   LEU B 226      -3.523 -18.584  10.065  1.00 30.00           O
ATOM     67  CB  LEU B 226      -4.587 -16.155   8.061  1.00 30.00           C
ATOM     68  CG  LEU B 226      -4.521 -16.312   6.542  1.00 30.00           C
ATOM     69  CD1 LEU B 226      -5.914 -16.258   5.943  1.00 30.00           C
ATOM     70  CD2 LEU B 226      -3.643 -15.232   5.923  1.00 30.00           C
ATOM     71  N   ILE B 227      -4.476 -16.749  11.064  1.00 30.00           N
ATOM     72  CA  ILE B 227      -3.690 -16.855  12.301  1.00 30.00           C
ATOM     73  C   ILE B 227      -3.952 -18.198  12.969  1.00 30.00           C
ATOM     74  O   ILE B 227      -3.000 -18.870  13.405  1.00 30.00           O
ATOM     75  CB  ILE B 227      -3.668 -15.644  13.279  1.00 30.00           C
ATOM     76  CG1 ILE B 227      -5.016 -15.359  13.918  1.00 30.00           C
ATOM     77  CG2 ILE B 227      -3.107 -14.410  12.588  1.00 30.00           C
ATOM     78  CD1 ILE B 227      -4.926 -14.502  15.160  1.00 30.00           C
ATOM     79  N   MET B 228      -5.227 -18.566  13.017  1.00 30.00           N
ATOM     80  CA  MET B 228      -5.663 -19.830  13.617  1.00 30.00           C
ATOM     81  C   MET B 228      -4.991 -21.002  12.916  1.00 30.00           C
ATOM     82  O   MET B 228      -4.472 -21.923  13.572  1.00 30.00           O
ATOM     83  CB  MET B 228      -7.192 -19.985  13.595  1.00 30.00           C
ATOM     84  CG  MET B 228      -7.960 -19.118  14.591  1.00 30.00           C
ATOM     85  SD  MET B 228      -7.511 -19.311  16.335  1.00 30.00           S
ATOM     86  CE  MET B 228      -7.660 -21.082  16.598  1.00 30.00           C
TER
"""
def exercise_02():
  """
  Generate ncs operators
  """
  from mmtbx.ncs.ncs import generate_ncs_ops
  assert generate_ncs_ops(symmetry='d7')[0].max_operators()==14
  assert generate_ncs_ops(symmetry='c7')[0].max_operators()==7
  assert generate_ncs_ops(symmetry='o')[0].max_operators()==24
  assert generate_ncs_ops(symmetry='t')[0].max_operators()==12
  assert generate_ncs_ops(symmetry='i')[0].max_operators()==60


def exercise_03():
  """
  Verify that there are no errors processing the write command
  No inception of the output is done. Just making sure it does not break
  """
  pdb_inp = pdb.input(source_info=None, lines=pdb_str_1)
  transform_info = pdb_inp.process_MTRIX_records()
  transforms_obj = iotbx.ncs.input(
    hierarchy=pdb_inp.construct_hierarchy())
  pdb_inp = pdb.input(source_info=None, lines=pdb_str_1)
  transforms_obj.get_ncs_info_as_spec()

def exercise_04():
  """Test MTRIX record processing"""
  expected  = """\
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
SCALE1      0.100000  0.000000  0.000000        0.00000
SCALE2      0.000000  0.100000  0.000000        0.00000
SCALE3      0.000000  0.000000  0.100000        0.00000
ATOM      1   N  ILE A  40       1.000   1.000   1.000  1.00162.33           C
ATOM      2  CA  LEU A  40      94.618  -5.253  91.582  1.00 87.10           C
TER
ATOM      3   C  ARG B  40      62.395  51.344  80.786  1.00107.25           C
TER
HETATM    4  C1  EDO A  40      39.954  51.526  72.372  0.33 60.93           C
ATOM      5   N  ILE C  40       1.000  -1.000   1.000  1.00162.33           C
ATOM      6  CA  LEU C  40      94.618 -91.582  -5.253  1.00 87.10           C
TER
ATOM      7   C  ARG D  40      62.395 -80.786  51.344  1.00107.25           C
TER
HETATM    8  C1  EDO C  40      39.954 -72.372  51.526  0.33 60.93           C
ATOM      9   N  ILE E  40       1.000   1.000  -1.000  1.00162.33           C
ATOM     10  CA  LEU E  40      91.582  -5.253 -94.618  1.00 87.10           C
TER
ATOM     11   C  ARG F  40      80.786  51.344 -62.395  1.00107.25           C
TER
HETATM   12  C1  EDO E  40      72.372  51.526 -39.954  0.33 60.93           C
ATOM     13   N  ILE G  40      -1.000   1.000   1.000  1.00162.33           C
ATOM     14  CA  LEU G  40       5.253  94.618  91.582  1.00 87.10           C
TER
ATOM     15   C  ARG H  40     -51.344  62.395  80.786  1.00107.25           C
TER
HETATM   16  C1  EDO G  40     -51.526  39.954  72.372  0.33 60.93           C
ATOM     17   N  ILE I  40       1.000   1.000  -1.000  1.00162.33           C
ATOM     18  CA  LEU I  40      91.582  -5.253 -94.618  1.00 87.10           C
TER
ATOM     19   C  ARG J  40      80.786  51.344 -62.395  1.00107.25           C
TER
HETATM   20  C1  EDO I  40      72.372  51.526 -39.954  0.33 60.93           C
ATOM     21   N  ILE K  40      -1.000   1.000  -1.000  1.00162.33           C
ATOM     22  CA  LEU K  40       5.253  91.582 -94.618  1.00 87.10           C
TER
ATOM     23   C  ARG L  40     -51.344  80.786 -62.395  1.00107.25           C
TER
HETATM   24  C1  EDO K  40     -51.526  72.372 -39.954  0.33 60.93           C
ATOM     25   N  ILE M  40      -0.366   1.366   1.000  1.00162.33           C
ATOM     26  CA  LEU M  40      51.858  79.315  91.582  1.00 87.10           C
TER
ATOM     27   C  ARG N  40     -13.268  79.708  80.786  1.00107.25           C
TER
HETATM   28  C1  EDO M  40     -24.646  60.364  72.372  0.33 60.93           C
ATOM     29   N  ILE O  40      -1.366   0.366   1.000  1.00162.33           C
ATOM     30  CA  LEU O  40     -42.760  84.568  91.582  1.00 87.10           C
TER
ATOM     31   C  ARG P  40     -75.663  28.364  80.786  1.00107.25           C
TER
HETATM   32  C1  EDO O  40     -64.600   8.838  72.372  0.33 60.93           C
ATOM     33   N  ILE Q  40       1.000   1.500   1.000  1.00162.33           C
ATOM     34  CA  LEU Q  40      94.618  -4.753  91.582  1.00 87.10           C
TER
ATOM     35   C  ARG R  40      62.395  51.844  80.786  1.00107.25           C
TER
HETATM   36  C1  EDO Q  40      39.954  52.026  72.372  0.33 60.93           C
ATOM     37   N  ILE S  40      -1.366   0.366   0.000  1.00162.33           C
ATOM     38  CA  LEU S  40     -42.760  84.568  90.582  1.00 87.10           C
TER
ATOM     39   C  ARG T  40     -75.663  28.364  79.786  1.00107.25           C
TER
HETATM   40  C1  EDO S  40     -64.600   8.838  71.372  0.33 60.93           C
END
"""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str_3)
  model = mmtbx.model.manager(pdb_inp, expand_with_mtrix=False)
  model.expand_with_BIOMT_records()
  assert not show_diff(expected, model.model_as_pdb())

def exercise_05():
  """Test MTRIX record processing"""
  cau_expected_results  = """\
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
SCALE1      0.100000  0.000000  0.000000        0.00000
SCALE2      0.000000  0.100000  0.000000        0.00000
SCALE3      0.000000  0.000000  0.100000        0.00000
ATOM      1   N  ILE A  40       1.000   1.000   1.000  1.00162.33           C
ATOM      2  CA  LEU A  40      94.618  -5.253  91.582  1.00 87.10           C
TER
ATOM      3   C  ARG B  40      62.395  51.344  80.786  1.00107.25           C
TER
HETATM    4  C1  EDO A  40      39.954  51.526  72.372  0.33 60.93           C
ATOM      5   N  ILE C  40       1.000  -1.000   1.000  1.00162.33           C
ATOM      6  CA  LEU C  40      94.618 -91.582  -5.253  1.00 87.10           C
TER
ATOM      7   C  ARG D  40      62.395 -80.786  51.344  1.00107.25           C
TER
HETATM    8  C1  EDO C  40      39.954 -72.372  51.526  0.33 60.93           C
ATOM      9   N  ILE E  40      -0.366   1.366   1.000  1.00162.33           C
ATOM     10  CA  LEU E  40      51.858  79.315  91.582  1.00 87.10           C
TER
ATOM     11   C  ARG F  40     -13.268  79.708  80.786  1.00107.25           C
TER
HETATM   12  C1  EDO E  40     -24.646  60.364  72.372  0.33 60.93           C
ATOM     13   N  ILE G  40      -1.366   0.366   1.000  1.00162.33           C
ATOM     14  CA  LEU G  40     -42.760  84.568  91.582  1.00 87.10           C
TER
ATOM     15   C  ARG H  40     -75.663  28.364  80.786  1.00107.25           C
TER
HETATM   16  C1  EDO G  40     -64.600   8.838  72.372  0.33 60.93           C
ATOM     17   N  ILE I  40       1.000   1.500   1.000  1.00162.33           C
ATOM     18  CA  LEU I  40      94.618  -4.753  91.582  1.00 87.10           C
TER
ATOM     19   C  ARG J  40      62.395  51.844  80.786  1.00107.25           C
TER
HETATM   20  C1  EDO I  40      39.954  52.026  72.372  0.33 60.93           C
END
"""
  # use MTRIX data
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str_3)
  model = mmtbx.model.manager(pdb_inp)
  assert not show_diff(cau_expected_results, model.model_as_pdb())

def exercise_06():
  """ Test that when building bio-molecule and then finding NCS relations
  from it, we get the same rotation and translation"""
  pdb_strings = [pdb_str_4, pdb_str_5]
  for method,pdb_string in enumerate(pdb_strings):
    pdb_inp = pdb.input(source_info=None, lines=pdb_string)
    model = mmtbx.model.manager(pdb_inp, expand_with_mtrix=False)
    crystal_symmetry = model.crystal_symmetry()
    # The exact transforms from pdb_string
    r1_expected = matrix.sqr(
      [0.309017, -0.951057, 0.0,0.951057, 0.309017,-0.0,0.0,0.0,1.0])
    r2_expected = matrix.sqr(
      [-0.809017,-0.587785,0.0,0.587785,-0.809017,-0.0,0.0,0.0,1.0])
    t1_expected = matrix.col([0,0,7])
    t2_expected = matrix.col([0,0,0])
    # Look at biomt records retrieved from PDB file
    if method == 0:
      rec = model._model_input.process_BIOMT_records()
      model.expand_with_BIOMT_records()
      h = model.get_hierarchy()
    else:
      rec = model._model_input.process_MTRIX_records()
      model.expand_with_MTRIX_records()
      h = model.get_hierarchy()
    r1 = rec.r[1]
    r2 = rec.r[2]
    t1 = rec.t[1]
    t2 = rec.t[2]
    assert approx_equal(r1, r1_expected, eps=0.001)
    assert approx_equal(t1, t1_expected, eps=0.1)
    assert approx_equal(r2, r2_expected, eps=0.001)
    assert approx_equal(t2, t2_expected, eps=0.1)
    # Look at the rotation and translation found by the NCS search
    s = h.as_pdb_string(crystal_symmetry=crystal_symmetry)
    ncs_obj = ncs.input(hierarchy=pdb.input(
      source_info=None, lines=s).construct_hierarchy())
    nrgl = ncs_obj.get_ncs_restraints_group_list()
    assert approx_equal(r1_expected, nrgl[0].copies[0].r, eps=0.001)
    assert approx_equal(t1_expected, nrgl[0].copies[0].t, eps=0.1)
    assert approx_equal(r2_expected, nrgl[0].copies[1].r, eps=0.001)
    assert approx_equal(t2_expected, nrgl[0].copies[1].t, eps=0.1)
    if method == 0:
      assert nrgl.get_n_groups() == 1
    elif method == 1:
      assert nrgl.get_n_groups() == 2

def exercise_08():
  """
  Test for MTRIX record when copies already present in file
  """
  pdb_inp = pdb.input(source_info=None, lines=pdb_str_8)
  model = mmtbx.model.manager(pdb_inp)
  assert model.get_number_of_atoms() == 7

if(__name__=='__main__'):
  exercise_02()
  exercise_03()
  exercise_04()
  exercise_05()
  exercise_06()
  exercise_08()
