from __future__ import absolute_import, division, print_function
from mmtbx.secondary_structure.build import ss_idealization as ssb
import iotbx.pdb
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from six.moves import zip


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
ATOM      3  C   ALA A   1      -0.114  -3.466  -2.684  1.00  0.00
ATOM      4  O   ALA A   1      -0.327  -3.827  -3.840  1.00  0.00
ATOM      5  CB  ALA A   1       1.308  -4.763  -1.090  1.00  0.00           C
ATOM      6  N   CYS A   2      -1.028  -2.938  -1.882  1.00  0.00
ATOM      7  CA  CYS A   2      -2.395  -2.743  -2.332  1.00  0.00
ATOM      8  C   CYS A   2      -2.396  -1.855  -3.579  1.00  0.00
ATOM      9  O   CYS A   2      -3.059  -2.167  -4.567  1.00  0.00
ATOM     10  CB  CYS A   2      -3.238  -2.101  -1.237  1.00  0.00           C
ATOM     11  SG  CYS A   2      -3.482  -3.154   0.212  1.00  0.00           S
ATOM     12  N   GLU A   3      -1.646  -0.767  -3.491  1.00  0.00
ATOM     13  CA  GLU A   3      -1.551   0.168  -4.599  1.00  0.00
ATOM     14  C   GLU A   3      -1.044  -0.568  -5.841  1.00  0.00
ATOM     15  O   GLU A   3      -1.601  -0.419  -6.927  1.00  0.00
ATOM     16  CB  GLU A   3      -0.610   1.317  -4.257  1.00  0.00           C
ATOM     17  CG  GLU A   3      -1.140   2.237  -3.169  1.00  0.00           C
ATOM     18  CD  GLU A   3      -0.179   3.362  -2.837  1.00  0.00           C
ATOM     19  OE1 GLU A   3       0.954   3.354  -3.365  1.00  0.00           O
ATOM     20  OE2 GLU A   3      -0.554   4.255  -2.048  1.00  0.00           O
ATOM     21  N   ASP A   4       0.008  -1.348  -5.639  1.00  0.00
ATOM     22  CA  ASP A   4       0.597  -2.109  -6.728  1.00  0.00
ATOM     23  C   ASP A   4      -0.466  -3.023  -7.340  1.00  0.00
ATOM     24  O   ASP A   4      -0.611  -3.085  -8.559  1.00  0.00
ATOM     25  CB  ASP A   4       1.773  -2.941  -6.231  1.00  0.00           C
ATOM     26  CG  ASP A   4       2.495  -3.656  -7.355  1.00  0.00           C
ATOM     27  OD1 ASP A   4       3.222  -2.985  -8.118  1.00  0.00           O
ATOM     28  OD2 ASP A   4       2.337  -4.889  -7.477  1.00  0.00           O
ATOM     29  N   GLY A   5      -1.184  -3.711  -6.463  1.00  0.00
ATOM     30  CA  GLY A   5      -2.231  -4.619  -6.901  1.00  0.00
ATOM     31  C   GLY A   5      -3.253  -3.847  -7.737  1.00  0.00
ATOM     32  O   GLY A   5      -3.647  -4.296  -8.813  1.00  0.00
ATOM     33  N   PHE A   6      -3.654  -2.699  -7.211  1.00  0.00
ATOM     34  CA  PHE A   6      -4.623  -1.860  -7.896  1.00  0.00
ATOM     35  C   PHE A   6      -4.090  -1.499  -9.284  1.00  0.00
ATOM     36  O   PHE A   6      -4.809  -1.602 -10.276  1.00  0.00
ATOM     37  CB  PHE A   6      -4.896  -0.590  -7.098  1.00  0.00           C
ATOM     38  CG  PHE A   6      -5.613  -0.838  -5.803  1.00  0.00           C
ATOM     39  CD1 PHE A   6      -6.985  -1.027  -5.781  1.00  0.00           C
ATOM     40  CD2 PHE A   6      -4.916  -0.882  -4.606  1.00  0.00           C
ATOM     41  CE1 PHE A   6      -7.649  -1.255  -4.589  1.00  0.00           C
ATOM     42  CE2 PHE A   6      -5.574  -1.110  -3.411  1.00  0.00           C
ATOM     43  CZ  PHE A   6      -6.942  -1.297  -3.404  1.00  0.00           C
ATOM     44  N   ILE A   7      -2.831  -1.084  -9.309  1.00  0.00
ATOM     45  CA  ILE A   7      -2.192  -0.708 -10.559  1.00  0.00
ATOM     46  C   ILE A   7      -2.243  -1.890 -11.529  1.00  0.00
ATOM     47  O   ILE A   7      -2.600  -1.727 -12.695  1.00  0.00
ATOM     48  CB  ILE A   7      -0.732  -0.284 -10.349  1.00  0.00           C
ATOM     49  CG1 ILE A   7       0.091  -1.437  -9.762  1.00  0.00           C
ATOM     50  CG2 ILE A   7      -0.673   0.933  -9.437  1.00  0.00           C
ATOM     51  CD1 ILE A   7       1.586  -1.180  -9.736  1.00  0.00           C
ATOM     52  N   HIS A   8      -1.881  -3.055 -11.012  1.00  0.00
ATOM     53  CA  HIS A   8      -1.882  -4.264 -11.817  1.00  0.00
ATOM     54  C   HIS A   8      -3.285  -4.496 -12.382  1.00  0.00
ATOM     55  O   HIS A   8      -3.442  -4.772 -13.570  1.00  0.00
ATOM     56  CB  HIS A   8      -1.451  -5.466 -10.986  1.00  0.00           C
ATOM     57  CG  HIS A   8      -0.031  -5.390 -10.518  1.00  0.00           C
ATOM     58  ND1 HIS A   8       1.041  -5.612 -11.354  1.00  0.00           N
ATOM     59  CD2 HIS A   8       0.493  -5.116  -9.299  1.00  0.00           C
ATOM     60  CE1 HIS A   8       2.164  -5.479 -10.671  1.00  0.00           C
ATOM     61  NE2 HIS A   8       1.859  -5.178  -9.422  1.00  0.00           N
ATOM     62  N   LYS A   9      -4.269  -4.376 -11.503  1.00  0.00
ATOM     63  CA  LYS A   9      -5.653  -4.568 -11.898  1.00  0.00
ATOM     64  C   LYS A   9      -6.000  -3.590 -13.022  1.00  0.00
ATOM     65  O   LYS A   9      -6.590  -3.978 -14.029  1.00  0.00
ATOM     66  CB  LYS A   9      -6.588  -4.352 -10.714  1.00  0.00           C
ATOM     67  CG  LYS A   9      -6.476  -5.423  -9.642  1.00  0.00           C
ATOM     68  CD  LYS A   9      -7.452  -5.175  -8.504  1.00  0.00           C
ATOM     69  CE  LYS A   9      -7.325  -6.235  -7.421  1.00  0.00           C
ATOM     70  NZ  LYS A   9      -8.267  -5.998  -6.293  1.00  0.00           N
ATOM     71  N   MET A  10      -5.617  -2.338 -12.812  1.00  0.00
ATOM     72  CA  MET A  10      -5.879  -1.301 -13.795  1.00  0.00
ATOM     73  C   MET A  10      -5.242  -1.695 -15.130  1.00  0.00
ATOM     74  O   MET A  10      -5.880  -1.605 -16.177  1.00  0.00
ATOM     75  CB  MET A  10      -5.315   0.036 -13.329  1.00  0.00           C
ATOM     76  CG  MET A  10      -6.072   0.644 -12.159  1.00  0.00           C
ATOM     77  SD  MET A  10      -5.738  -0.196 -10.599  1.00  0.00           S
ATOM     78  CE  MET A  10      -6.442   0.965  -9.431  1.00  0.00           C
ATOM     79  N   LEU A  11      -3.991  -2.124 -15.047  1.00  0.00
ATOM     80  CA  LEU A  11      -3.261  -2.533 -16.235  1.00  0.00
ATOM     81  C   LEU A  11      -4.024  -3.658 -16.937  1.00  0.00
ATOM     82  O   LEU A  11      -4.213  -3.623 -18.151  1.00  0.00
ATOM     83  CB  LEU A  11      -1.859  -3.009 -15.869  1.00  0.00           C
ATOM     84  CG  LEU A  11      -0.901  -1.945 -15.320  1.00  0.00           C
ATOM     85  CD1 LEU A  11       0.404  -2.587 -14.875  1.00  0.00           C
ATOM     86  CD2 LEU A  11      -0.638  -0.851 -16.346  1.00  0.00           C
ATOM     87  N   ASN A  12      -4.443  -4.631 -16.140  1.00  0.00
ATOM     88  CA  ASN A  12      -5.183  -5.765 -16.669  1.00  0.00
ATOM     89  C   ASN A  12      -6.440  -5.262 -17.382  1.00  0.00
ATOM     90  O   ASN A  12      -6.738  -5.685 -18.497  1.00  0.00
ATOM     91  CB  ASN A  12      -5.566  -6.727 -15.557  1.00  0.00           C
ATOM     92  CG  ASN A  12      -4.368  -7.481 -15.014  1.00  0.00           C
ATOM     93  OD1 ASN A  12      -3.501  -7.919 -15.770  1.00  0.00           O
ATOM     94  ND2 ASN A  12      -4.314  -7.637 -13.696  1.00  0.00           N
ATOM     95  N   GLN A  13      -7.144  -4.364 -16.707  1.00  0.00
ATOM     96  CA  GLN A  13      -8.362  -3.798 -17.261  1.00  0.00
ATOM     97  C   GLN A  13      -8.047  -3.133 -18.603  1.00  0.00
ATOM     98  O   GLN A  13      -8.755  -3.341 -19.586  1.00  0.00
ATOM     99  CB  GLN A  13      -8.966  -2.779 -16.302  1.00  0.00           C
ATOM    100  CG  GLN A  13     -10.288  -2.188 -16.771  1.00  0.00           C
ATOM    101  CD  GLN A  13     -10.889  -1.230 -15.763  1.00  0.00           C
ATOM    102  OE1 GLN A  13     -10.303  -0.967 -14.712  1.00  0.00           O
ATOM    103  NE2 GLN A  13     -12.066  -0.701 -16.077  1.00  0.00           N
ATOM    104  N   PRO A  14      -6.981  -2.345 -18.600  1.00  0.00
ATOM    105  CA  PRO A  14      -6.562  -1.648 -19.804  1.00  0.00
ATOM    106  C   PRO A  14      -6.302  -2.668 -20.915  1.00  0.00
ATOM    107  O   PRO A  14      -6.758  -2.492 -22.043  1.00  0.00
ATOM    108  CB  PRO A  14      -5.265  -0.968 -19.395  1.00  0.00           C
ATOM    109  CG  PRO A  14      -5.393  -0.780 -17.926  1.00  0.00           C
ATOM    110  CD  PRO A  14      -6.119  -1.992 -17.421  1.00  0.00           C
ATOM    111  N   SER A  15      -5.570  -3.712 -20.555  1.00  0.00
ATOM    112  CA  SER A  15      -5.244  -4.761 -21.507  1.00  0.00
ATOM    113  C   SER A  15      -6.538  -5.354 -22.070  1.00  0.00
ATOM    114  O   SER A  15      -6.671  -5.527 -23.280  1.00  0.00
ATOM    115  CB  SER A  15      -4.417  -5.855 -20.843  1.00  0.00           C
ATOM    116  OG  SER A  15      -4.136  -6.904 -21.753  1.00  0.00           O
ATOM    117  N   ARG A  16      -7.458  -5.648 -21.163  1.00  0.00
ATOM    118  CA  ARG A  16      -8.737  -6.218 -21.553  1.00  0.00
ATOM    119  C   ARG A  16      -9.432  -5.275 -22.538  1.00  0.00
ATOM    120  O   ARG A  16      -9.931  -5.711 -23.574  1.00  0.00
ATOM    121  CB  ARG A  16      -9.625  -6.440 -20.333  1.00 10.00           C
ATOM    122  CG  ARG A  16      -9.062  -7.437 -19.324  1.00 10.00           C
ATOM    123  CD  ARG A  16      -8.986  -8.857 -19.879  1.00 10.00           C
ATOM    124  NE  ARG A  16     -10.308  -9.411 -20.184  1.00 10.00           N
ATOM    125  CZ  ARG A  16     -10.555 -10.473 -20.955  1.00 10.00           C
ATOM    126  NH1 ARG A  16      -9.582 -11.166 -21.548  1.00 10.00           N
ATOM    127  NH2 ARG A  16     -11.812 -10.855 -21.138  1.00 10.00           N
ATOM    128  N   THR A  17      -9.442  -4.000 -22.179  1.00  0.00
ATOM    129  CA  THR A  17     -10.067  -2.991 -23.017  1.00  0.00
ATOM    130  C   THR A  17      -9.418  -3.011 -24.403  1.00  0.00
ATOM    131  O   THR A  17     -10.112  -3.007 -25.418  1.00  0.00
ATOM    132  CB  THR A  17      -9.922  -1.586 -22.421  1.00  0.00           C
ATOM    133  OG1 THR A  17      -8.535  -1.249 -22.296  1.00  0.00           O
ATOM    134  CG2 THR A  17     -10.584  -1.520 -21.053  1.00  0.00           C
ATOM    135  N   TRP A  18      -8.093  -3.034 -24.400  1.00  0.00
ATOM    136  CA  TRP A  18      -7.342  -3.056 -25.643  1.00  0.00
ATOM    137  C   TRP A  18      -7.761  -4.276 -26.467  1.00  0.00
ATOM    138  O   TRP A  18      -8.021  -4.163 -27.663  1.00  0.00
ATOM    139  CB  TRP A  18      -5.844  -3.104 -25.369  1.00  0.00           C
ATOM    140  CG  TRP A  18      -5.314  -1.843 -24.760  1.00  0.00           C
ATOM    141  CD1 TRP A  18      -5.088  -1.602 -23.436  1.00  0.00           C
ATOM    142  CD2 TRP A  18      -4.942  -0.645 -25.454  1.00  0.00           C
ATOM    143  NE1 TRP A  18      -4.599  -0.330 -23.263  1.00  0.00           N
ATOM    144  CE2 TRP A  18      -4.500   0.278 -24.486  1.00  0.00           C
ATOM    145  CE3 TRP A  18      -4.941  -0.264 -26.800  1.00  0.00           C
ATOM    146  CZ2 TRP A  18      -4.061   1.557 -24.820  1.00  0.00           C
ATOM    147  CZ3 TRP A  18      -4.504   1.006 -27.130  1.00  0.00           C
ATOM    148  CH2 TRP A  18      -4.070   1.902 -26.144  1.00  0.00           C
ATOM    149  N   VAL A  19      -7.813  -5.415 -25.792  1.00  0.00
ATOM    150  CA  VAL A  19      -8.197  -6.655 -26.445  1.00  0.00
ATOM    151  C   VAL A  19      -9.587  -6.493 -27.064  1.00  0.00
ATOM    152  O   VAL A  19      -9.803  -6.858 -28.218  1.00  0.00
ATOM    153  CB  VAL A  19      -8.224  -7.833 -25.460  1.00  0.00           C
ATOM    154  CG1 VAL A  19      -8.745  -9.091 -26.144  1.00  0.00           C
ATOM    155  CG2 VAL A  19      -6.836  -8.068 -24.880  1.00  0.00           C
ATOM    156  N   TYR A  20     -10.493  -5.943 -26.268  1.00  0.00
ATOM    157  CA  TYR A  20     -11.856  -5.727 -26.723  1.00  0.00
ATOM    158  C   TYR A  20     -11.839  -4.848 -27.975  1.00  0.00
ATOM    159  O   TYR A  20     -12.504  -5.155 -28.963  1.00  0.00
ATOM    160  CB  TYR A  20     -12.690  -5.063 -25.634  1.00  0.00           C
ATOM    161  CG  TYR A  20     -12.961  -5.962 -24.448  1.00  0.00           C
ATOM    162  CD1 TYR A  20     -12.741  -7.332 -24.525  1.00  0.00           C
ATOM    163  CD2 TYR A  20     -13.437  -5.441 -23.252  1.00  0.00           C
ATOM    164  CE1 TYR A  20     -12.987  -8.157 -23.444  1.00  0.00           C
ATOM    165  CE2 TYR A  20     -13.686  -6.259 -22.165  1.00  0.00           C
ATOM    166  CZ  TYR A  20     -13.459  -7.615 -22.267  1.00  0.00           C
ATOM    167  OH  TYR A  20     -13.705  -8.433 -21.187  1.00  0.00           O
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
ATOM      1  N   ALA A   1     -34.967   0.361   0.480  1.00 30.00           N
ATOM      2  CA  ALA A   1     -33.700   0.830  -0.121  1.00 30.00           C
ATOM      3  C   ALA A   1     -32.491   0.227   0.567  1.00 30.00           C
ATOM      4  O   ALA A   1     -32.528   0.205   1.784  1.00 30.00           O
ATOM         CB  ALA A   1     -33.642   2.352  -0.067  1.00  0.00           C
ATOM      1  N   ALA A   2     -31.523  -0.190  -0.301  1.00 30.00           N
ATOM      2  CA  ALA A   2     -30.230  -0.801   0.073  1.00 30.00           C
ATOM      3  C   ALA A   2     -29.061   0.113  -0.241  1.00 30.00           C
ATOM      4  O   ALA A   2     -29.096   0.668  -1.324  1.00 30.00           O
ATOM         CB  ALA A   2     -30.071  -2.137  -0.644  1.00  0.00           C
ATOM      1  N   ALA A   3     -28.126   0.160   0.752  1.00 30.00           N
ATOM      2  CA  ALA A   3     -26.874   0.947   0.728  1.00 30.00           C
ATOM      3  C   ALA A   3     -25.649   0.058   0.644  1.00 30.00           C
ATOM      4  O   ALA A   3     -25.650  -0.919   1.370  1.00 30.00           O
ATOM         CB  ALA A   3     -26.807   1.836   1.964  1.00  0.00           C
ATOM      1  N   ALA A   4     -24.710   0.510  -0.238  1.00 30.00           N
ATOM      2  CA  ALA A   4     -23.409  -0.133  -0.523  1.00 30.00           C
ATOM      3  C   ALA A   4     -22.245   0.697  -0.017  1.00 30.00           C
ATOM      4  O   ALA A   4     -22.309   1.892  -0.239  1.00 30.00           O
ATOM         CB  ALA A   4     -23.279  -0.379  -2.021  1.00  0.00           C
ATOM      1  N   ALA A   5     -21.280  -0.042   0.605  1.00 30.00           N
ATOM      2  CA  ALA A   5     -20.026   0.484   1.186  1.00 30.00           C
ATOM      3  C   ALA A   5     -18.809   0.032   0.403  1.00 30.00           C
ATOM      4  O   ALA A   5     -18.794  -1.140   0.074  1.00 30.00           O
ATOM         CB  ALA A   5     -19.917   0.050   2.642  1.00  0.00           C
ATOM      1  N   ALA A   6     -17.893   1.024   0.201  1.00 30.00           N
ATOM      2  CA  ALA A   6     -16.605   0.883  -0.511  1.00 30.00           C
ATOM      3  C   ALA A   6     -15.423   1.013   0.428  1.00 30.00           C
ATOM      4  O   ALA A   6     -15.486   1.918   1.241  1.00 30.00           O
ATOM         CB  ALA A   6     -16.521   1.921  -1.624  1.00  0.00           C
ATOM      1  N   ALA A   7     -14.444   0.089   0.202  1.00 30.00           N
ATOM      2  CA  ALA A   7     -13.170  -0.024   0.945  1.00 30.00           C
ATOM      3  C   ALA A   7     -11.979   0.343   0.082  1.00 30.00           C
ATOM      4  O   ALA A   7     -11.980  -0.112  -1.048  1.00 30.00           O
ATOM         CB  ALA A   7     -13.020  -1.440   1.488  1.00  0.00           C
ATOM      1  N   ALA A   8     -11.064   1.128   0.723  1.00 30.00           N
ATOM      2  CA  ALA A   8      -9.799   1.631   0.147  1.00 30.00           C
ATOM      3  C   ALA A   8      -8.588   0.991   0.797  1.00 30.00           C
ATOM      4  O   ALA A   8      -8.623   0.900   2.011  1.00 30.00           O
ATOM         CB  ALA A   8      -9.741   3.148   0.287  1.00  0.00           C
ATOM      1  N   ALA A   9      -7.622   0.624  -0.094  1.00 30.00           N
ATOM      2  CA  ALA A   9      -6.327  -0.006   0.241  1.00 30.00           C
ATOM      3  C   ALA A   9      -5.160   0.925  -0.022  1.00 30.00           C
ATOM      4  O   ALA A   9      -5.197   1.540  -1.072  1.00 30.00           O
ATOM         CB  ALA A   9      -6.169  -1.299  -0.550  1.00  0.00           C
ATOM      1  N   ALA A  10      -4.222   0.916   0.970  1.00 30.00           N
ATOM      2  CA  ALA A  10      -2.972   1.704   0.988  1.00 30.00           C
ATOM      3  C   ALA A  10      -1.746   0.822   0.852  1.00 30.00           C
ATOM      4  O   ALA A  10      -1.745  -0.195   1.522  1.00 30.00           O
ATOM         CB  ALA A  10      -2.902   2.521   2.273  1.00  0.00           C
ATOM      1  N   ALA A  11      -0.808   1.323  -0.005  1.00 30.00           N
ATOM      2  CA  ALA A  11       0.492   0.698  -0.328  1.00 30.00           C
ATOM      3  C   ALA A  11       1.656   1.499   0.222  1.00 30.00           C
ATOM      4  O   ALA A  11       1.591   2.705   0.068  1.00 30.00           O
ATOM         CB  ALA A  11       0.620   0.538  -1.838  1.00  0.00           C
ATOM      1  N   ALA A  12       2.623   0.727   0.800  1.00 30.00           N
ATOM      2  CA  ALA A  12       3.878   1.220   1.406  1.00 30.00           C
ATOM      3  C   ALA A  12       5.094   0.813   0.597  1.00 30.00           C
ATOM      4  O   ALA A  12       5.109  -0.338   0.202  1.00 30.00           O
ATOM         CB  ALA A  12       3.990   0.704   2.836  1.00  0.00           C
ATOM      1  N   ALA A  13       6.009   1.816   0.450  1.00 30.00           N
ATOM      2  CA  ALA A  13       7.295   1.716  -0.272  1.00 30.00           C
ATOM      3  C   ALA A  13       8.480   1.793   0.672  1.00 30.00           C
ATOM      4  O   ALA A  13       8.417   2.651   1.534  1.00 30.00           O
ATOM         CB  ALA A  13       7.377   2.815  -1.324  1.00  0.00           C
ATOM      1  N   ALA A  14       9.458   0.884   0.392  1.00 30.00           N
ATOM      2  CA  ALA A  14      10.734   0.730   1.124  1.00 30.00           C
ATOM      3  C   ALA A  14      11.923   1.146   0.282  1.00 30.00           C
ATOM      4  O   ALA A  14      11.920   0.756  -0.872  1.00 30.00           O
ATOM         CB  ALA A  14      10.885  -0.714   1.587  1.00  0.00           C
ATOM      1  N   ALA A  15      12.839   1.894   0.964  1.00 30.00           N
ATOM      2  CA  ALA A  15      14.103   2.430   0.416  1.00 30.00           C
ATOM      3  C   ALA A  15      15.315   1.754   1.027  1.00 30.00           C
ATOM      4  O   ALA A  15      15.283   1.594   2.233  1.00 30.00           O
ATOM         CB  ALA A  15      14.160   3.936   0.642  1.00  0.00           C
ATOM      1  N   ALA A  16      16.280   1.439   0.114  1.00 30.00           N
ATOM      2  CA  ALA A  16      17.575   0.791   0.411  1.00 30.00           C
ATOM      3  C   ALA A  16      18.742   1.737   0.199  1.00 30.00           C
ATOM      4  O   ALA A  16      18.702   2.410  -0.815  1.00 30.00           O
ATOM         CB  ALA A  16      17.733  -0.454  -0.453  1.00  0.00           C
ATOM      1  N   ALA A  17      19.681   1.672   1.187  1.00 30.00           N
ATOM      2  CA  ALA A  17      20.932   2.458   1.247  1.00 30.00           C
ATOM      3  C   ALA A  17      22.157   1.587   1.058  1.00 30.00           C
ATOM      4  O   ALA A  17      22.160   0.533   1.670  1.00 30.00           O
ATOM         CB  ALA A  17      21.003   3.201   2.576  1.00  0.00           C
ATOM      1  N   ALA A  18      23.093   2.136   0.230  1.00 30.00           N
ATOM      2  CA  ALA A  18      24.393   1.531  -0.131  1.00 30.00           C
ATOM      3  C   ALA A  18      25.558   2.300   0.462  1.00 30.00           C
ATOM      4  O   ALA A  18      25.492   3.513   0.377  1.00 30.00           O
ATOM         CB  ALA A  18      24.518   1.457  -1.648  1.00  0.00           C
ATOM      1  N   ALA A  19      26.526   1.497   0.993  1.00 30.00           N
ATOM      2  CA  ALA A  19      27.782   1.955   1.625  1.00 30.00           C
ATOM      3  C   ALA A  19      28.997   1.596   0.792  1.00 30.00           C
ATOM      4  O   ALA A  19      29.011   0.469   0.331  1.00 30.00           O
ATOM         CB  ALA A  19      27.897   1.359   3.022  1.00  0.00           C
ATOM      1  N   ALA A  20      29.911   2.606   0.699  1.00 30.00           N
ATOM      2  CA  ALA A  20      31.196   2.548  -0.029  1.00 30.00           C
ATOM      3  C   ALA A  20      32.382   2.572   0.915  1.00 30.00           C
ATOM      4  O   ALA A  20      32.321   3.379   1.825  1.00 30.00           O
ATOM         CB  ALA A  20      31.275   3.705  -1.017  1.00  0.00           C
ATOM      1  N   ALA A  21      33.361   1.681   0.582  1.00 30.00           N
ATOM      2  CA  ALA A  21      34.637   1.486   1.303  1.00 30.00           C
ATOM      3  C   ALA A  21      35.825   1.950   0.482  1.00 30.00           C
ATOM      4  O   ALA A  21      35.820   1.626  -0.691  1.00 30.00           O
ATOM         CB  ALA A  21      34.791   0.018   1.682  1.00  0.00           C
ATOM      1  N   ALA A  22      36.741   2.659   1.205  1.00 30.00           N
ATOM      2  CA  ALA A  22      38.004   3.225   0.685  1.00 30.00           C
ATOM      3  C   ALA A  22      39.218   2.517   1.254  1.00 30.00           C
ATOM      4  O   ALA A  22      39.188   2.289   2.450  1.00 30.00           O
ATOM         CB  ALA A  22      38.061   4.716   0.996  1.00  0.00           C
TER
"""

def exercise_00(prefix="exercise_00"):
  "Build poly-ALA alpha helix"
  ph = ssb.secondary_structure_from_sequence(
    alpha_helix_template, "".join(["A"]*12))
  ph.atoms().reset_i_seq()
  #ph.write_pdb_file(file_name="%s_result.pdb"%prefix)
  sites_1 = ph.atoms().extract_xyz()
  sites_2 = iotbx.pdb.input(source_info=None,
    lines=alpha_helix_answer).construct_hierarchy().atoms().extract_xyz()
  #of = open("1.pdb","w")
  #print >> of, alpha_helix_answer
  #of.close()
  assert sites_1.size()==sites_2.size(), [sites_1.size(),sites_2.size()]
  d1 = flex.sqrt((sites_1 - sites_2).dot())
  rmsd = ssb.calculate_rmsd_smart(ph, iotbx.pdb.input(source_info=None,
    lines=alpha_helix_answer).construct_hierarchy())
  assert rmsd < 0.1, rmsd

def exercise_01(prefix="exercise_01"):
  "Build poly-ALA beta strand"
  ph = ssb.secondary_structure_from_sequence(
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
  except Exception as e:
    assert str(e) == "pdb_str should contain at least 2 residues"

def exercise_secondary_structure_from_sequence():
  pdb_inp = iotbx.pdb.input(source_info=None, lines=correct_answer)
  cs = pdb_inp.xray_structure_simple().crystal_symmetry()
  correct_h = pdb_inp.construct_hierarchy()
  test_h = ssb.secondary_structure_from_sequence(alpha_pdb_str,
    sequence="ACEDGFIHKMLNQPSRTWVY")
  test_h.atoms().reset_serial()
  # correct_h.write_pdb_file("correct_h.pdb")
  # test_h.write_pdb_file("test_h.pdb")
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
  assert test_h.as_str() == correct_h.as_str()
  #correct_h.write_pdb_file(file_name="tst_1_ex_last_correct.pdb")
  #test_h.write_pdb_file(file_name="tst_1_ex_last_out.pdb")
  # ...and this fails ?
  # assert correct_h.is_similar_hierarchy(other=test_h)

  assert approx_equal(test_h.atoms().extract_xyz(),
                      correct_h.atoms().extract_xyz(), eps=0.002)
  try: ssb.secondary_structure_from_sequence(alpha_pdb_str, sequence="")
  except Exception as e:
    assert str(e) == "sequence should contain at least one residue."


def exercise():
  exercise_00()
  exercise_01()
  exercise_r_t_matrices()
  exercise_secondary_structure_from_sequence()

if (__name__ == "__main__"):
  exercise()
