from __future__ import division
from mmtbx.secondary_structure import build as ssb
import iotbx.pdb
from libtbx.test_utils import approx_equal
from scitbx.math import superpose
from scitbx.array_family import flex


t_pdb_str = """\
ATOM      1  N   ALA     2       1.643  -2.366  -1.408  1.00
ATOM      3  CA  ALA     2       1.280  -3.608  -2.069  1.00
ATOM      6  CB  ALA     2       1.361  -4.762  -1.068  1.00
ATOM     10  C   ALA     2      -0.114  -3.466  -2.684  1.00
ATOM     11  O   ALA     2      -0.327  -3.827  -3.840  1.00
"""

alpha_pdb_str = """\
ATOM      1  N   ALA     2       1.643  -2.366  -1.408  1.00
ATOM      3  CA  ALA     2       1.280  -3.608  -2.069  1.00
ATOM      6  CB  ALA     2       1.361  -4.762  -1.068  1.00
ATOM     10  C   ALA     2      -0.114  -3.466  -2.684  1.00
ATOM     11  O   ALA     2      -0.327  -3.827  -3.840  1.00
ATOM     12  N   ALA     3      -1.028  -2.938  -1.882  1.00
ATOM     14  CA  ALA     3      -2.395  -2.743  -2.333  1.00
ATOM     17  CB  ALA     3      -3.228  -2.150  -1.194  1.00
ATOM     21  C   ALA     3      -2.396  -1.855  -3.579  1.00
ATOM     22  O   ALA     3      -3.059  -2.167  -4.567  1.00
"""

correct_answer = """\
ATOM      1  N   ALA A   1       1.643  -2.366  -1.408  1.00  0.00
ATOM      2  CA  ALA A   1       1.280  -3.608  -2.069  1.00  0.00
ATOM      3  CB  ALA A   1       1.361  -4.762  -1.068  1.00  0.00
ATOM      4  C   ALA A   1      -0.114  -3.466  -2.684  1.00  0.00
ATOM      5  O   ALA A   1      -0.327  -3.827  -3.840  1.00  0.00
ATOM      6  N   CYS A   2      -1.028  -2.938  -1.882  1.00  0.00
ATOM      7  CA  CYS A   2      -2.395  -2.743  -2.332  1.00  0.00
ATOM      8  CB  CYS A   2      -3.228  -2.150  -1.194  1.00  0.00
ATOM      9  C   CYS A   2      -2.396  -1.855  -3.579  1.00  0.00
ATOM     10  O   CYS A   2      -3.059  -2.167  -4.567  1.00  0.00
ATOM     11  SG  CYS A   2      -3.433  -3.241   0.233  1.00  0.00           S
ATOM     12  N   GLU A   3      -1.646  -0.767  -3.491  1.00  0.00
ATOM     13  CA  GLU A   3      -1.551   0.168  -4.599  1.00  0.00
ATOM     14  CB  GLU A   3      -0.646   1.337  -4.205  1.00  0.00
ATOM     15  C   GLU A   3      -1.044  -0.568  -5.841  1.00  0.00
ATOM     16  O   GLU A   3      -1.601  -0.419  -6.927  1.00  0.00
ATOM     17  CG  GLU A   3      -1.217   2.231  -3.116  1.00  0.00           C
ATOM     18  CD  GLU A   3      -0.297   3.383  -2.759  1.00  0.00           C
ATOM     19  OE1 GLU A   3       0.848   3.121  -2.333  1.00  0.00           O
ATOM     20  OE2 GLU A   3      -0.719   4.550  -2.904  1.00  0.00           O
ATOM     21  N   ASP A   4       0.008  -1.348  -5.639  1.00  0.00
ATOM     22  CA  ASP A   4       0.597  -2.109  -6.728  1.00  0.00
ATOM     23  CB  ASP A   4       1.808  -2.887  -6.211  1.00  0.00
ATOM     24  C   ASP A   4      -0.466  -3.023  -7.340  1.00  0.00
ATOM     25  O   ASP A   4      -0.611  -3.085  -8.559  1.00  0.00
ATOM     26  CG  ASP A   4       2.490  -3.691  -7.310  1.00  0.00           C
ATOM     27  OD1 ASP A   4       3.032  -3.072  -8.251  1.00  0.00           O
ATOM     28  OD2 ASP A   4       2.483  -4.938  -7.231  1.00  0.00           O
ATOM     29  N   GLY A   5      -1.184  -3.711  -6.463  1.00  0.00
ATOM     30  CA  GLY A   5      -2.231  -4.619  -6.901  1.00  0.00
ATOM     31  C   GLY A   5      -3.253  -3.847  -7.737  1.00  0.00
ATOM     32  O   GLY A   5      -3.647  -4.296  -8.813  1.00  0.00
ATOM     33  N   PHE A   6      -3.654  -2.699  -7.211  1.00  0.00
ATOM     34  CA  PHE A   6      -4.623  -1.860  -7.896  1.00  0.00
ATOM     35  CB  PHE A   6      -4.919  -0.623  -7.045  1.00  0.00
ATOM     36  C   PHE A   6      -4.090  -1.499  -9.284  1.00  0.00
ATOM     37  O   PHE A   6      -4.809  -1.602 -10.276  1.00  0.00
ATOM     38  CG  PHE A   6      -5.579  -0.924  -5.731  1.00  0.00           C
ATOM     39  CD1 PHE A   6      -6.947  -1.123  -5.654  1.00  0.00           C
ATOM     40  CD2 PHE A   6      -4.827  -1.013  -4.571  1.00  0.00           C
ATOM     41  CE1 PHE A   6      -7.554  -1.405  -4.443  1.00  0.00           C
ATOM     42  CE2 PHE A   6      -5.427  -1.294  -3.358  1.00  0.00           C
ATOM     43  CZ  PHE A   6      -6.792  -1.491  -3.295  1.00  0.00           C
ATOM     44  N   ILE A   7      -2.831  -1.084  -9.309  1.00  0.00
ATOM     45  CA  ILE A   7      -2.192  -0.708 -10.559  1.00  0.00
ATOM     46  CB  ILE A   7      -0.761  -0.243 -10.281  1.00  0.00
ATOM     47  C   ILE A   7      -2.243  -1.890 -11.529  1.00  0.00
ATOM     48  O   ILE A   7      -2.600  -1.727 -12.695  1.00  0.00
ATOM     49  CG1 ILE A   7       0.086  -1.363  -9.667  1.00  0.00           C
ATOM     50  CG2 ILE A   7      -0.745   0.991  -9.393  1.00  0.00           C
ATOM     51  CD1 ILE A   7       1.571  -1.063  -9.579  1.00  0.00           C
ATOM     52  N   HIS A   8      -1.881  -3.055 -11.012  1.00  0.00
ATOM     53  CA  HIS A   8      -1.882  -4.264 -11.817  1.00  0.00
ATOM     54  CB  HIS A   8      -1.391  -5.441 -10.972  1.00  0.00
ATOM     55  C   HIS A   8      -3.285  -4.496 -12.382  1.00  0.00
ATOM     56  O   HIS A   8      -3.442  -4.772 -13.570  1.00  0.00
ATOM     57  CG  HIS A   8      -0.003  -5.277 -10.440  1.00  0.00           C
ATOM     58  ND1 HIS A   8       1.112  -5.285 -11.250  1.00  0.00           N
ATOM     59  CD2 HIS A   8       0.450  -5.092  -9.178  1.00  0.00           C
ATOM     60  CE1 HIS A   8       2.192  -5.114 -10.509  1.00  0.00           C
ATOM     61  NE2 HIS A   8       1.818  -4.994  -9.248  1.00  0.00           N
ATOM     62  N   LYS A   9      -4.269  -4.376 -11.503  1.00  0.00
ATOM     63  CA  LYS A   9      -5.653  -4.568 -11.898  1.00  0.00
ATOM     64  CB  LYS A   9      -6.561  -4.400 -10.678  1.00  0.00
ATOM     65  C   LYS A   9      -6.000  -3.590 -13.022  1.00  0.00
ATOM     66  O   LYS A   9      -6.590  -3.978 -14.029  1.00  0.00
ATOM     67  CG  LYS A   9      -6.405  -5.485  -9.625  1.00  0.00           C
ATOM     68  CD  LYS A   9      -7.268  -5.210  -8.405  1.00  0.00           C
ATOM     69  CE  LYS A   9      -7.123  -6.307  -7.362  1.00  0.00           C
ATOM     70  NZ  LYS A   9      -7.907  -6.016  -6.132  1.00  0.00           N
ATOM     71  N   MET A  10      -5.617  -2.338 -12.812  1.00  0.00
ATOM     72  CA  MET A  10      -5.879  -1.301 -13.795  1.00  0.00
ATOM     73  CB  MET A  10      -5.358   0.040 -13.274  1.00  0.00
ATOM     74  C   MET A  10      -5.242  -1.695 -15.130  1.00  0.00
ATOM     75  O   MET A  10      -5.880  -1.605 -16.177  1.00  0.00
ATOM     76  CG  MET A  10      -6.112   0.579 -12.070  1.00  0.00           C
ATOM     77  SD  MET A  10      -5.421   2.120 -11.440  1.00  0.00           S
ATOM     78  CE  MET A  10      -6.482   2.416 -10.027  1.00  0.00           C
ATOM     79  N   LEU A  11      -3.991  -2.124 -15.047  1.00  0.00
ATOM     80  CA  LEU A  11      -3.261  -2.533 -16.235  1.00  0.00
ATOM     81  CB  LEU A  11      -1.841  -2.948 -15.844  1.00  0.00
ATOM     82  C   LEU A  11      -4.024  -3.658 -16.937  1.00  0.00
ATOM     83  O   LEU A  11      -4.213  -3.623 -18.151  1.00  0.00
ATOM     84  CG  LEU A  11      -0.925  -1.853 -15.287  1.00  0.00           C
ATOM     85  CD1 LEU A  11       0.395  -2.447 -14.820  1.00  0.00           C
ATOM     86  CD2 LEU A  11      -0.685  -0.751 -16.311  1.00  0.00           C
ATOM     87  N   ASN A  12      -4.443  -4.631 -16.140  1.00  0.00
ATOM     88  CA  ASN A  12      -5.183  -5.765 -16.669  1.00  0.00
ATOM     89  CB  ASN A  12      -5.506  -6.738 -15.534  1.00  0.00
ATOM     90  C   ASN A  12      -6.440  -5.262 -17.382  1.00  0.00
ATOM     91  O   ASN A  12      -6.738  -5.685 -18.497  1.00  0.00
ATOM     92  CG  ASN A  12      -4.285  -7.482 -15.031  1.00  0.00           C
ATOM     93  OD1 ASN A  12      -3.288  -7.612 -15.741  1.00  0.00           O
ATOM     94  ND2 ASN A  12      -4.355  -7.973 -13.799  1.00  0.00           N
ATOM     95  N   GLN A  13      -7.144  -4.364 -16.707  1.00  0.00
ATOM     96  CA  GLN A  13      -8.362  -3.798 -17.261  1.00  0.00
ATOM     97  CB  GLN A  13      -8.975  -2.821 -16.256  1.00  0.00
ATOM     98  C   GLN A  13      -8.047  -3.133 -18.603  1.00  0.00
ATOM     99  O   GLN A  13      -8.755  -3.341 -19.586  1.00  0.00
ATOM    100  CG  GLN A  13     -10.338  -2.274 -16.676  1.00  0.00           C
ATOM    101  CD  GLN A  13     -10.995  -1.444 -15.590  1.00  0.00           C
ATOM    102  OE1 GLN A  13     -10.482  -1.342 -14.476  1.00  0.00           O
ATOM    103  NE2 GLN A  13     -12.137  -0.847 -15.910  1.00  0.00           N
ATOM    104  N   PRO A  14      -6.981  -2.345 -18.600  1.00  0.00
ATOM    105  CA  PRO A  14      -6.562  -1.648 -19.804  1.00  0.00
ATOM    106  CB  PRO A  14      -5.330  -0.794 -19.497  1.00  0.00
ATOM    107  C   PRO A  14      -6.302  -2.668 -20.915  1.00  0.00
ATOM    108  O   PRO A  14      -6.758  -2.492 -22.043  1.00  0.00
ATOM    109  CG  PRO A  14      -5.536  -0.534 -18.013  1.00  0.00           C
ATOM    110  CD  PRO A  14      -6.105  -1.788 -17.417  1.00  0.00           C
ATOM    111  N   SER A  15      -5.570  -3.712 -20.555  1.00  0.00
ATOM    112  CA  SER A  15      -5.244  -4.761 -21.507  1.00  0.00
ATOM    113  CB  SER A  15      -4.367  -5.814 -20.827  1.00  0.00
ATOM    114  C   SER A  15      -6.538  -5.354 -22.070  1.00  0.00
ATOM    115  O   SER A  15      -6.671  -5.527 -23.280  1.00  0.00
ATOM    116  OG  SER A  15      -4.069  -6.880 -21.724  1.00  0.00           O
ATOM    117  N   ARG A  16      -7.458  -5.648 -21.163  1.00  0.00
ATOM    118  CA  ARG A  16      -8.737  -6.218 -21.553  1.00  0.00
ATOM    119  CB  ARG A  16      -9.580  -6.482 -20.305  1.00  0.00
ATOM    120  C   ARG A  16      -9.432  -5.275 -22.538  1.00  0.00
ATOM    121  O   ARG A  16      -9.931  -5.711 -23.574  1.00  0.00
ATOM    122  CG  ARG A  16      -9.007  -7.529 -19.354  1.00 10.00           C
ATOM    123  CD  ARG A  16      -9.141  -8.953 -19.888  1.00 10.00           C
ATOM    124  NE  ARG A  16     -10.536  -9.396 -19.958  1.00 10.00           N
ATOM    125  CZ  ARG A  16     -10.974 -10.526 -20.519  1.00 10.00           C
ATOM    126  NH1 ARG A  16     -10.150 -11.401 -21.096  1.00 10.00           N
ATOM    127  NH2 ARG A  16     -12.274 -10.789 -20.503  1.00 10.00           N
ATOM    128  N   THR A  17      -9.442  -4.000 -22.179  1.00  0.00
ATOM    129  CA  THR A  17     -10.067  -2.991 -23.017  1.00  0.00
ATOM    130  CB  THR A  17      -9.955  -1.623 -22.341  1.00  0.00
ATOM    131  C   THR A  17      -9.418  -3.011 -24.403  1.00  0.00
ATOM    132  O   THR A  17     -10.112  -3.007 -25.418  1.00  0.00
ATOM    133  OG1 THR A  17      -8.584  -1.270 -22.134  1.00  0.00           O
ATOM    134  CG2 THR A  17     -10.697  -1.590 -21.014  1.00  0.00           C
ATOM    135  N   TRP A  18      -8.093  -3.034 -24.400  1.00  0.00
ATOM    136  CA  TRP A  18      -7.342  -3.056 -25.643  1.00  0.00
ATOM    137  CB  TRP A  18      -5.844  -3.047 -25.335  1.00  0.00
ATOM    138  C   TRP A  18      -7.761  -4.276 -26.467  1.00  0.00
ATOM    139  O   TRP A  18      -8.021  -4.163 -27.663  1.00  0.00
ATOM    140  CG  TRP A  18      -5.350  -1.753 -24.766  1.00  0.00           C
ATOM    141  CD1 TRP A  18      -5.192  -1.450 -23.445  1.00  0.00           C
ATOM    142  CD2 TRP A  18      -4.955  -0.584 -25.495  1.00  0.00           C
ATOM    143  NE1 TRP A  18      -4.724  -0.166 -23.306  1.00  0.00           N
ATOM    144  CE2 TRP A  18      -4.570   0.387 -24.550  1.00  0.00           C
ATOM    145  CE3 TRP A  18      -4.890  -0.265 -26.856  1.00  0.00           C
ATOM    146  CZ2 TRP A  18      -4.127   1.655 -24.921  1.00  0.00           C
ATOM    147  CZ3 TRP A  18      -4.450   0.994 -27.222  1.00  0.00           C
ATOM    148  CH2 TRP A  18      -4.073   1.938 -26.258  1.00  0.00           C
ATOM    149  N   VAL A  19      -7.813  -5.415 -25.792  1.00  0.00
ATOM    150  CA  VAL A  19      -8.197  -6.655 -26.445  1.00  0.00
ATOM    151  CB  VAL A  19      -8.138  -7.803 -25.436  1.00  0.00
ATOM    152  C   VAL A  19      -9.587  -6.493 -27.064  1.00  0.00
ATOM    153  O   VAL A  19      -9.803  -6.858 -28.218  1.00  0.00
ATOM    154  CG1 VAL A  19      -8.661  -9.102 -26.086  1.00  0.00           C
ATOM    155  CG2 VAL A  19      -6.724  -8.021 -24.917  1.00  0.00           C
ATOM    156  N   TYR A  20     -10.493  -5.943 -26.268  1.00  0.00
ATOM    157  CA  TYR A  20     -11.856  -5.727 -26.723  1.00  0.00
ATOM    158  CB  TYR A  20     -12.680  -5.113 -25.590  1.00  0.00
ATOM    159  C   TYR A  20     -11.839  -4.848 -27.975  1.00  0.00
ATOM    160  O   TYR A  20     -12.504  -5.155 -28.963  1.00  0.00
ATOM    161  CG  TYR A  20     -12.847  -6.011 -24.387  1.00  0.00           C
ATOM    162  CD1 TYR A  20     -13.096  -7.370 -24.538  1.00  0.00           C
ATOM    163  CD2 TYR A  20     -12.745  -5.504 -23.098  1.00  0.00           C
ATOM    164  CE1 TYR A  20     -13.241  -8.197 -23.440  1.00  0.00           C
ATOM    165  CE2 TYR A  20     -12.889  -6.323 -21.993  1.00  0.00           C
ATOM    166  CZ  TYR A  20     -13.136  -7.668 -22.171  1.00  0.00           C
ATOM    167  OH  TYR A  20     -13.279  -8.488 -21.074  1.00  0.00           O
TER
"""

alpha_helix_template="""\
ATOM      1  N   ALA     1      -1.204  -0.514   0.643  1.00  0.00           N
ATOM      2  CA  ALA     1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  CB  ALA     1       0.808   0.860   0.974  1.00  0.00           C
ATOM      4  C   ALA     1       0.866  -1.134  -0.537  1.00  0.00           C
ATOM      5  O   ALA     1       1.432  -1.035  -1.625  1.00  0.00           O
ATOM      6  N   ALA     2       0.965  -2.213   0.234  1.00  0.00           N
ATOM      7  CA  ALA     2       1.761  -3.368  -0.162  1.00  0.00           C
ATOM      8  CB  ALA     2       1.751  -4.434   0.936  1.00  0.00           C
ATOM      9  C   ALA     2       1.258  -3.963  -1.472  1.00  0.00           C
ATOM     10  O   ALA     2       2.047  -4.366  -2.327  1.00  0.00           O
"""

alpha_helix_answer="""\
ATOM      1  N   ALA     1      -1.204  -0.514   0.643  1.00  0.00           N
ATOM      2  CA  ALA     1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  CB  ALA     1       0.808   0.860   0.974  1.00  0.00           C
ATOM      4  C   ALA     1       0.866  -1.134  -0.537  1.00  0.00           C
ATOM      5  O   ALA     1       1.432  -1.035  -1.625  1.00  0.00           O
ATOM      6  N   ALA     2       0.965  -2.213   0.234  1.00  0.00           N
ATOM      7  CA  ALA     2       1.761  -3.368  -0.162  1.00  0.00           C
ATOM      8  CB  ALA     2       1.751  -4.434   0.936  1.00  0.00           C
ATOM      9  C   ALA     2       1.258  -3.963  -1.472  1.00  0.00           C
ATOM     10  O   ALA     2       2.047  -4.366  -2.327  1.00  0.00           O
ATOM     11  N   ALA     3      -0.062  -4.016  -1.625  1.00  0.00           N
ATOM     12  CA  ALA     3      -0.673  -4.562  -2.830  1.00  0.00           C
ATOM     13  CB  ALA     3      -2.199  -4.560  -2.711  1.00  0.00           C
ATOM     14  C   ALA     3      -0.245  -3.782  -4.068  1.00  0.00           C
ATOM     15  O   ALA     3       0.011  -4.364  -5.123  1.00  0.00           O
ATOM     16  N   ALA     4      -0.168  -2.462  -3.934  1.00  0.00           N
ATOM     17  CA  ALA     4       0.229  -1.600  -5.040  1.00  0.00           C
ATOM     18  CB  ALA     4       0.169  -0.128  -4.627  1.00  0.00           C
ATOM     19  C   ALA     4       1.628  -1.946  -5.537  1.00  0.00           C
ATOM     20  O   ALA     4       1.888  -1.949  -6.740  1.00  0.00           O
ATOM     21  N   ALA     5       2.529  -2.239  -4.602  1.00  0.00           N
ATOM     22  CA  ALA     5       3.902  -2.587  -4.943  1.00  0.00           C
ATOM     23  CB  ALA     5       4.732  -2.811  -3.677  1.00  0.00           C
ATOM     24  C   ALA     5       3.954  -3.827  -5.827  1.00  0.00           C
ATOM     25  O   ALA     5       4.744  -3.899  -6.767  1.00  0.00           O
ATOM     26  N   ALA     6       3.104  -4.804  -5.520  1.00  0.00           N
ATOM     27  CA  ALA     6       3.052  -6.043  -6.286  1.00  0.00           C
ATOM     28  CB  ALA     6       2.036  -7.012  -5.678  1.00  0.00           C
ATOM     29  C   ALA     6       2.707  -5.774  -7.747  1.00  0.00           C
ATOM     30  O   ALA     6       3.266  -6.394  -8.653  1.00  0.00           O
ATOM     31  N   ALA     7       1.783  -4.845  -7.972  1.00  0.00           N
ATOM     32  CA  ALA     7       1.361  -4.493  -9.322  1.00  0.00           C
ATOM     33  CB  ALA     7       0.245  -3.447  -9.284  1.00  0.00           C
ATOM     34  C   ALA     7       2.533  -3.976 -10.149  1.00  0.00           C
ATOM     35  O   ALA     7       2.658  -4.294 -11.332  1.00  0.00           O
ATOM     36  N   ALA     8       3.390  -3.178  -9.521  1.00  0.00           N
ATOM     37  CA  ALA     8       4.552  -2.615 -10.197  1.00  0.00           C
ATOM     38  CB  ALA     8       5.326  -1.686  -9.260  1.00  0.00           C
ATOM     39  C   ALA     8       5.472  -3.714 -10.720  1.00  0.00           C
ATOM     40  O   ALA     8       6.014  -3.613 -11.821  1.00  0.00           O
ATOM     41  N   ALA     9       5.645  -4.766  -9.923  1.00  0.00           N
ATOM     42  CA  ALA     9       6.497  -5.884 -10.304  1.00  0.00           C
ATOM     43  CB  ALA     9       6.565  -6.919  -9.179  1.00  0.00           C
ATOM     44  C   ALA     9       6.006  -6.541 -11.590  1.00  0.00           C
ATOM     45  O   ALA     9       6.801  -6.921 -12.448  1.00  0.00           O
ATOM     46  N   ALA    10       4.688  -6.669 -11.716  1.00  0.00           N
ATOM     47  CA  ALA    10       4.088  -7.281 -12.896  1.00  0.00           C
ATOM     48  CB  ALA    10       2.567  -7.360 -12.750  1.00  0.00           C
ATOM     49  C   ALA    10       4.452  -6.512 -14.162  1.00  0.00           C
ATOM     50  O   ALA    10       4.723  -7.107 -15.206  1.00  0.00           O
ATOM     51  N   ALA    11       4.458  -5.187 -14.064  1.00  0.00           N
ATOM     52  CA  ALA    11       4.789  -4.335 -15.199  1.00  0.00           C
ATOM     53  CB  ALA    11       4.655  -2.857 -14.824  1.00  0.00           C
ATOM     54  C   ALA    11       6.198  -4.617 -15.709  1.00  0.00           C
ATOM     55  O   ALA    11       6.441  -4.639 -16.916  1.00  0.00           O
ATOM     56  N   ALA    12       7.126  -4.834 -14.781  1.00  0.00           N
ATOM     57  CA  ALA    12       8.512  -5.115 -15.134  1.00  0.00           C
ATOM     58  CB  ALA    12       9.371  -5.259 -13.875  1.00  0.00           C
ATOM     59  C   ALA    12       8.622  -6.375 -15.987  1.00  0.00           C
ATOM     60  O   ALA    12       9.403  -6.427 -16.937  1.00  0.00           O
TER
END
"""

beta_strand_template="""\
ATOM      1  N   ALA A   1     -34.967   0.361   0.480  1.00 30.00           N
ATOM      2  CA  ALA A   1     -33.700   0.830  -0.121  1.00 30.00           C
ATOM      3  C   ALA A   1     -32.491   0.227   0.567  1.00 30.00           C
ATOM      4  O   ALA A   1     -32.528   0.205   1.784  1.00 30.00           O
ATOM      5  CB  ALA A   1     -33.586   2.364  -0.139  1.00 30.00           C
ATOM      6  N   ALA A   2     -31.529  -0.207  -0.287  1.00 30.00           N
ATOM      7  CA  ALA A   2     -30.265  -0.830   0.121  1.00 30.00           C
ATOM      8  C   ALA A   2     -29.092   0.104  -0.268  1.00 30.00           C
ATOM      9  O   ALA A   2     -29.024   0.723  -1.359  1.00 30.00           O
ATOM     10  CB  ALA A   2     -29.837  -2.178  -0.470  1.00 30.00           C
"""

beta_strand_answer="""\
ATOM      1  N   ALA     1     -34.967   0.361   0.480  1.00 30.00           N
ATOM      2  CA  ALA     1     -33.700   0.830  -0.121  1.00 30.00           C
ATOM      3  C   ALA     1     -32.491   0.227   0.567  1.00 30.00           C
ATOM      4  O   ALA     1     -32.528   0.205   1.784  1.00 30.00           O
ATOM      5  CB  ALA     1     -33.586   2.364  -0.139  1.00 30.00           C
ATOM      6  N   ALA     2     -31.523  -0.190  -0.301  1.00 30.00           N
ATOM      7  CA  ALA     2     -30.230  -0.801   0.073  1.00 30.00           C
ATOM      8  C   ALA     2     -29.061   0.113  -0.241  1.00 30.00           C
ATOM      9  O   ALA     2     -29.096   0.668  -1.324  1.00 30.00           O
ATOM     10  CB  ALA     2     -30.014  -2.176  -0.583  1.00 30.00           C
ATOM     11  N   ALA     3     -28.126   0.160   0.752  1.00 30.00           N
ATOM     12  CA  ALA     3     -26.874   0.947   0.728  1.00 30.00           C
ATOM     13  C   ALA     3     -25.649   0.058   0.644  1.00 30.00           C
ATOM     14  O   ALA     3     -25.650  -0.919   1.370  1.00 30.00           O
ATOM     15  CB  ALA     3     -26.753   1.901   1.928  1.00 30.00           C
ATOM     16  N   ALA     4     -24.710   0.510  -0.238  1.00 30.00           N
ATOM     17  CA  ALA     4     -23.409  -0.133  -0.523  1.00 30.00           C
ATOM     18  C   ALA     4     -22.245   0.697  -0.017  1.00 30.00           C
ATOM     19  O   ALA     4     -22.309   1.892  -0.239  1.00 30.00           O
ATOM     20  CB  ALA     4     -23.220  -0.450  -2.017  1.00 30.00           C
ATOM     21  N   ALA     5     -21.280  -0.042   0.605  1.00 30.00           N
ATOM     22  CA  ALA     5     -20.026   0.484   1.186  1.00 30.00           C
ATOM     23  C   ALA     5     -18.809   0.032   0.403  1.00 30.00           C
ATOM     24  O   ALA     5     -18.794  -1.140   0.074  1.00 30.00           O
ATOM     25  CB  ALA     5     -19.864   0.120   2.671  1.00 30.00           C
ATOM     26  N   ALA     6     -17.893   1.024   0.201  1.00 30.00           N
ATOM     27  CA  ALA     6     -16.605   0.883  -0.511  1.00 30.00           C
ATOM     28  C   ALA     6     -15.423   1.013   0.428  1.00 30.00           C
ATOM     29  O   ALA     6     -15.486   1.918   1.241  1.00 30.00           O
ATOM     30  CB  ALA     6     -16.463   1.875  -1.679  1.00 30.00           C
ATOM     31  N   ALA     7     -14.444   0.089   0.202  1.00 30.00           N
ATOM     32  CA  ALA     7     -13.170  -0.024   0.945  1.00 30.00           C
ATOM     33  C   ALA     7     -11.979   0.343   0.082  1.00 30.00           C
ATOM     34  O   ALA     7     -11.980  -0.112  -1.048  1.00 30.00           O
ATOM     35  CB  ALA     7     -12.965  -1.419   1.560  1.00 30.00           C
ATOM     36  N   ALA     8     -11.064   1.128   0.723  1.00 30.00           N
ATOM     37  CA  ALA     8      -9.799   1.631   0.147  1.00 30.00           C
ATOM     38  C   ALA     8      -8.588   0.991   0.797  1.00 30.00           C
ATOM     39  O   ALA     8      -8.623   0.900   2.011  1.00 30.00           O
ATOM     40  CB  ALA     8      -9.685   3.164   0.216  1.00 30.00           C
ATOM     41  N   ALA     9      -7.622   0.624  -0.094  1.00 30.00           N
ATOM     42  CA  ALA     9      -6.327  -0.006   0.241  1.00 30.00           C
ATOM     43  C   ALA     9      -5.160   0.925  -0.022  1.00 30.00           C
ATOM     44  O   ALA     9      -5.197   1.540  -1.072  1.00 30.00           O
ATOM     45  CB  ALA     9      -6.112  -1.342  -0.492  1.00 30.00           C
ATOM     46  N   ALA    10      -4.222   0.916   0.970  1.00 30.00           N
ATOM     47  CA  ALA    10      -2.972   1.704   0.988  1.00 30.00           C
ATOM     48  C   ALA    10      -1.746   0.822   0.852  1.00 30.00           C
ATOM     49  O   ALA    10      -1.745  -0.195   1.522  1.00 30.00           O
ATOM     50  CB  ALA    10      -2.848   2.588   2.241  1.00 30.00           C
ATOM     51  N   ALA    11      -0.808   1.323  -0.005  1.00 30.00           N
ATOM     52  CA  ALA    11       0.492   0.698  -0.328  1.00 30.00           C
ATOM     53  C   ALA    11       1.656   1.499   0.222  1.00 30.00           C
ATOM     54  O   ALA    11       1.591   2.705   0.068  1.00 30.00           O
ATOM     55  CB  ALA    11       0.678   0.466  -1.837  1.00 30.00           C
ATOM     56  N   ALA    12       2.623   0.727   0.800  1.00 30.00           N
ATOM     57  CA  ALA    12       3.878   1.220   1.406  1.00 30.00           C
ATOM     58  C   ALA    12       5.094   0.813   0.597  1.00 30.00           C
ATOM     59  O   ALA    12       5.109  -0.338   0.202  1.00 30.00           O
ATOM     60  CB  ALA    12       4.043   0.772   2.869  1.00 30.00           C
ATOM     61  N   ALA    13       6.009   1.816   0.450  1.00 30.00           N
ATOM     62  CA  ALA    13       7.295   1.716  -0.272  1.00 30.00           C
ATOM     63  C   ALA    13       8.480   1.793   0.672  1.00 30.00           C
ATOM     64  O   ALA    13       8.417   2.651   1.534  1.00 30.00           O
ATOM     65  CB  ALA    13       7.435   2.772  -1.381  1.00 30.00           C
ATOM     66  N   ALA    14       9.458   0.884   0.392  1.00 30.00           N
ATOM     67  CA  ALA    14      10.734   0.730   1.124  1.00 30.00           C
ATOM     68  C   ALA    14      11.923   1.146   0.282  1.00 30.00           C
ATOM     69  O   ALA    14      11.920   0.756  -0.872  1.00 30.00           O
ATOM     70  CB  ALA    14      10.940  -0.698   1.659  1.00 30.00           C
ATOM     71  N   ALA    15      12.839   1.894   0.964  1.00 30.00           N
ATOM     72  CA  ALA    15      14.103   2.430   0.416  1.00 30.00           C
ATOM     73  C   ALA    15      15.315   1.754   1.027  1.00 30.00           C
ATOM     74  O   ALA    15      15.283   1.594   2.233  1.00 30.00           O
ATOM     75  CB  ALA    15      14.216   3.956   0.571  1.00 30.00           C
ATOM     76  N   ALA    16      16.280   1.439   0.114  1.00 30.00           N
ATOM     77  CA  ALA    16      17.575   0.791   0.411  1.00 30.00           C
ATOM     78  C   ALA    16      18.742   1.737   0.199  1.00 30.00           C
ATOM     79  O   ALA    16      18.702   2.410  -0.815  1.00 30.00           O
ATOM     80  CB  ALA    16      17.790  -0.500  -0.397  1.00 30.00           C
ATOM     81  N   ALA    17      19.681   1.672   1.187  1.00 30.00           N
ATOM     82  CA  ALA    17      20.932   2.458   1.247  1.00 30.00           C
ATOM     83  C   ALA    17      22.157   1.587   1.058  1.00 30.00           C
ATOM     84  O   ALA    17      22.160   0.533   1.670  1.00 30.00           O
ATOM     85  CB  ALA    17      21.057   3.271   2.548  1.00 30.00           C
ATOM     86  N   ALA    18      23.093   2.136   0.230  1.00 30.00           N
ATOM     87  CA  ALA    18      24.393   1.531  -0.131  1.00 30.00           C
ATOM     88  C   ALA    18      25.558   2.300   0.462  1.00 30.00           C
ATOM     89  O   ALA    18      25.492   3.513   0.377  1.00 30.00           O
ATOM     90  CB  ALA    18      24.576   1.385  -1.651  1.00 30.00           C
ATOM     91  N   ALA    19      26.526   1.497   0.993  1.00 30.00           N
ATOM     92  CA  ALA    19      27.782   1.955   1.625  1.00 30.00           C
ATOM     93  C   ALA    19      28.997   1.596   0.792  1.00 30.00           C
ATOM     94  O   ALA    19      29.011   0.469   0.331  1.00 30.00           O
ATOM     95  CB  ALA    19      27.950   1.425   3.059  1.00 30.00           C
ATOM     96  N   ALA    20      29.911   2.606   0.699  1.00 30.00           N
ATOM     97  CA  ALA    20      31.196   2.548  -0.029  1.00 30.00           C
ATOM     98  C   ALA    20      32.382   2.572   0.915  1.00 30.00           C
ATOM     99  O   ALA    20      32.321   3.379   1.825  1.00 30.00           O
ATOM    100  CB  ALA    20      31.333   3.665  -1.077  1.00 30.00           C
ATOM    101  N   ALA    21      33.361   1.681   0.582  1.00 30.00           N
ATOM    102  CA  ALA    21      34.637   1.486   1.303  1.00 30.00           C
ATOM    103  C   ALA    21      35.825   1.950   0.482  1.00 30.00           C
ATOM    104  O   ALA    21      35.820   1.626  -0.691  1.00 30.00           O
ATOM    105  CB  ALA    21      34.845   0.031   1.755  1.00 30.00           C
ATOM    106  N   ALA    22      36.741   2.659   1.205  1.00 30.00           N
ATOM    107  CA  ALA    22      38.004   3.225   0.685  1.00 30.00           C
ATOM    108  C   ALA    22      39.218   2.517   1.254  1.00 30.00           C
ATOM    109  O   ALA    22      39.188   2.289   2.450  1.00 30.00           O
ATOM    110  CB  ALA    22      38.117   4.740   0.926  1.00 30.00           C
TER
"""

def exercise_00(prefix="exercise_00"):
  "Build poly-ALA alpha helix"
  ph = ssb.make_ss_structure_from_sequence(
    alpha_helix_template, "".join(["A"]*12))
  #ph.write_pdb_file(file_name="%s_result.pdb"%prefix)
  sites_1 = ph.atoms().extract_xyz()
  sites_2 = iotbx.pdb.input(source_info=None,
    lines=alpha_helix_answer).construct_hierarchy().atoms().extract_xyz()
  assert sites_1.size()==sites_2.size()
  d1 = flex.sqrt((sites_1 - sites_2).dot())
  assert d1.min_max_mean().as_tuple()[2]<0.02
  lsq_fit_obj = superpose.least_squares_fit(reference_sites = sites_1,
                                            other_sites     = sites_2)
  moving_sites = lsq_fit_obj.other_sites_best_fit()
  d2 = flex.sqrt((sites_1 - moving_sites).dot())
  assert d2.min_max_mean().as_tuple()[2]<0.005

def exercise_01(prefix="exercise_01"):
  "Build poly-ALA beta strand"
  ph = ssb.make_ss_structure_from_sequence(
    beta_strand_template, "".join(["A"]*22))
  sites_1 = ph.atoms().extract_xyz()
  sites_2 = iotbx.pdb.input(source_info=None,
    lines=beta_strand_answer).construct_hierarchy().atoms().extract_xyz()
  assert sites_1.size()==sites_2.size()
  assert approx_equal(sites_1, sites_2, eps=0.002)

def exercise_r_t_matrices():
  r,t = ssb.get_r_t_matrices_from_structure(alpha_pdb_str)
  assert approx_equal(r.elems,
                      (-0.02358, 0.96052,  0.27720,
                       -0.86374, -0.15919, 0.47811,
                        0.50337, -0.22815, 0.83340),
                      eps = 0.0001)
  assert approx_equal(t.elems,
                      (1.6740, -1.2223, -2.0756),
                      eps = 0.0001)
  try: ssb.get_r_t_matrices_from_structure(t_pdb_str)
  except Exception, e:
    assert str(e) == "pdb_str should contain at least 2 residues"

def exercise_ss_structure_from_sequence():
  pdb_inp = iotbx.pdb.input(source_info=None, lines=correct_answer)
  cs = pdb_inp.xray_structure_simple().crystal_symmetry()
  correct_h = pdb_inp.construct_hierarchy()
  test_h = ssb.make_ss_structure_from_sequence(alpha_pdb_str,
    sequence="ACEDGFIHKMLNQPSRTWVY")
  test_h.atoms().reset_serial()
  # why this passes...
  for m1, m2 in zip(correct_h.models(), test_h.models()):
    m1.is_similar_hierarchy(  other=m2)
    m1.is_identical_hierarchy(other=m2)
    for c1, c2 in zip(m1.chains(), m2.chains()):
      c1.is_similar_hierarchy(  other=c2)
      c1.is_identical_hierarchy(other=c2)
      for rg1, rg2 in zip(c1.residue_groups(), c2.residue_groups()):
        rg1.is_similar_hierarchy(  other=rg2)
        rg1.is_identical_hierarchy(other=rg2)
  f1 = correct_h.extract_xray_structure(crystal_symmetry=cs).structure_factors(
    algorithm="direct", d_min=2).f_calc().data()
  f2 = test_h.extract_xray_structure(crystal_symmetry=cs).structure_factors(
    algorithm="direct", d_min=2).f_calc().data()
  assert approx_equal(f1, f2)
  assert test_h.as_str() == correct_h.as_str()
  #correct_h.write_pdb_file(file_name="tst_1_ex_last_correct.pdb")
  #test_h.write_pdb_file(file_name="tst_1_ex_last_out.pdb")
  # ...and this fails ?
  #assert correct_h.is_similar_hierarchy(other=test_h)

  assert approx_equal(test_h.atoms().extract_xyz(),
                      correct_h.atoms().extract_xyz(), eps=0.002)
  try: ssb.make_ss_structure_from_sequence(alpha_pdb_str, sequence="")
  except Exception, e:
    assert str(e) == "sequence should contain at least one residue."


def exercise():
  exercise_00()
  exercise_01()
  exercise_r_t_matrices()
  exercise_ss_structure_from_sequence()

if (__name__ == "__main__"):
  exercise()
