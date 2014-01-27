from mmtbx.secondary_structure import build as ssb
import iotbx.pdb
from libtbx.test_utils import Exception_expected, approx_equal
from libtbx.utils import Sorry
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
ATOM      1  N   ALA     1       1.643  -2.366  -1.408  1.00  0.00
ATOM      2  CA  ALA     1       1.280  -3.608  -2.069  1.00  0.00
ATOM      3  CB  ALA     1       1.361  -4.762  -1.068  1.00  0.00
ATOM      4  C   ALA     1      -0.114  -3.466  -2.684  1.00  0.00
ATOM      5  O   ALA     1      -0.327  -3.827  -3.840  1.00  0.00
ATOM      6  N   CYS     2      -1.028  -2.938  -1.882  1.00  0.00
ATOM      7  CA  CYS     2      -2.395  -2.743  -2.332  1.00  0.00
ATOM      8  CB  CYS     2      -3.228  -2.150  -1.194  1.00  0.00
ATOM      9  C   CYS     2      -2.396  -1.855  -3.579  1.00  0.00
ATOM     10  O   CYS     2      -3.059  -2.167  -4.567  1.00  0.00
ATOM     11  SG  CYS     2      -3.315  -3.180   0.289  1.00  0.00           S
ATOM     12  N   GLU     3      -1.646  -0.767  -3.491  1.00  0.00
ATOM     13  CA  GLU     3      -1.551   0.168  -4.599  1.00  0.00
ATOM     14  CB  GLU     3      -0.646   1.337  -4.205  1.00  0.00
ATOM     15  C   GLU     3      -1.044  -0.568  -5.841  1.00  0.00
ATOM     16  O   GLU     3      -1.601  -0.419  -6.927  1.00  0.00
ATOM     17  CG  GLU     3      -1.142   2.140  -3.014  1.00  0.00           C
ATOM     18  CD  GLU     3      -0.228   3.298  -2.664  1.00  0.00           C
ATOM     19  OE1 GLU     3      -0.535   4.029  -1.699  1.00  0.00           O
ATOM     20  OE2 GLU     3       0.797   3.480  -3.355  1.00  0.00           O
ATOM     21  N   ASP     4       0.008  -1.348  -5.639  1.00  0.00
ATOM     22  CA  ASP     4       0.597  -2.109  -6.728  1.00  0.00
ATOM     23  CB  ASP     4       1.808  -2.887  -6.211  1.00  0.00
ATOM     24  C   ASP     4      -0.466  -3.023  -7.340  1.00  0.00
ATOM     25  O   ASP     4      -0.611  -3.085  -8.559  1.00  0.00
ATOM     26  CG  ASP     4       2.475  -3.712  -7.305  1.00  0.00           C
ATOM     27  OD1 ASP     4       3.387  -3.183  -7.976  1.00  0.00           O
ATOM     28  OD2 ASP     4       2.085  -4.884  -7.490  1.00  0.00           O
ATOM     29  N   GLY     5      -1.184  -3.711  -6.463  1.00  0.00
ATOM     30  CA  GLY     5      -2.231  -4.619  -6.901  1.00  0.00
ATOM     31  C   GLY     5      -3.253  -3.847  -7.737  1.00  0.00
ATOM     32  O   GLY     5      -3.647  -4.296  -8.813  1.00  0.00
ATOM     33  N   PHE     6      -3.654  -2.699  -7.211  1.00  0.00
ATOM     34  CA  PHE     6      -4.623  -1.860  -7.896  1.00  0.00
ATOM     35  CB  PHE     6      -4.919  -0.623  -7.045  1.00  0.00
ATOM     36  C   PHE     6      -4.090  -1.499  -9.284  1.00  0.00
ATOM     37  O   PHE     6      -4.809  -1.602 -10.276  1.00  0.00
ATOM     38  CG  PHE     6      -5.512  -0.929  -5.701  1.00  0.00           C
ATOM     39  CD1 PHE     6      -6.880  -1.080  -5.546  1.00  0.00           C
ATOM     40  CD2 PHE     6      -4.697  -1.068  -4.589  1.00  0.00           C
ATOM     41  CE1 PHE     6      -7.424  -1.365  -4.308  1.00  0.00           C
ATOM     42  CE2 PHE     6      -5.236  -1.354  -3.348  1.00  0.00           C
ATOM     43  CZ  PHE     6      -6.601  -1.503  -3.209  1.00  0.00           C
ATOM     44  N   ILE     7      -2.831  -1.084  -9.309  1.00  0.00
ATOM     45  CA  ILE     7      -2.192  -0.708 -10.559  1.00  0.00
ATOM     46  CB  ILE     7      -0.761  -0.243 -10.281  1.00  0.00
ATOM     47  C   ILE     7      -2.243  -1.890 -11.529  1.00  0.00
ATOM     48  O   ILE     7      -2.600  -1.727 -12.695  1.00  0.00
ATOM     49  CG1 ILE     7       0.038  -1.299  -9.514  1.00  0.00           C
ATOM     50  CG2 ILE     7      -0.750   1.085  -9.537  1.00  0.00           C
ATOM     51  CD1 ILE     7       1.515  -0.989  -9.354  1.00  0.00           C
ATOM     52  N   HIS     8      -1.881  -3.055 -11.012  1.00  0.00
ATOM     53  CA  HIS     8      -1.882  -4.264 -11.817  1.00  0.00
ATOM     54  CB  HIS     8      -1.391  -5.441 -10.972  1.00  0.00
ATOM     55  C   HIS     8      -3.285  -4.496 -12.382  1.00  0.00
ATOM     56  O   HIS     8      -3.442  -4.772 -13.570  1.00  0.00
ATOM     57  CG  HIS     8      -0.003  -5.277 -10.441  1.00  0.00           C
ATOM     58  ND1 HIS     8       1.116  -5.691 -11.131  1.00  0.00           N
ATOM     59  CD2 HIS     8       0.448  -4.736  -9.284  1.00  0.00           C
ATOM     60  CE1 HIS     8       2.195  -5.413 -10.422  1.00  0.00           C
ATOM     61  NE2 HIS     8       1.818  -4.834  -9.297  1.00  0.00           N
ATOM     62  N   LYS     9      -4.269  -4.376 -11.503  1.00  0.00
ATOM     63  CA  LYS     9      -5.653  -4.568 -11.898  1.00  0.00
ATOM     64  CB  LYS     9      -6.561  -4.400 -10.678  1.00  0.00
ATOM     65  C   LYS     9      -6.000  -3.590 -13.022  1.00  0.00
ATOM     66  O   LYS     9      -6.590  -3.978 -14.029  1.00  0.00
ATOM     67  CG  LYS     9      -6.297  -5.390  -9.557  1.00  0.00           C
ATOM     68  CD  LYS     9      -7.242  -5.178  -8.386  1.00  0.00           C
ATOM     69  CE  LYS     9      -6.971  -6.168  -7.265  1.00  0.00           C
ATOM     70  NZ  LYS     9      -7.889  -5.971  -6.111  1.00  0.00           N
ATOM     71  N   MET    10      -5.617  -2.338 -12.812  1.00  0.00
ATOM     72  CA  MET    10      -5.879  -1.301 -13.795  1.00  0.00
ATOM     73  CB  MET    10      -5.358   0.040 -13.274  1.00  0.00
ATOM     74  C   MET    10      -5.242  -1.695 -15.130  1.00  0.00
ATOM     75  O   MET    10      -5.880  -1.605 -16.177  1.00  0.00
ATOM     76  CG  MET    10      -6.005   0.503 -11.980  1.00  0.00           C
ATOM     77  SD  MET    10      -5.732  -0.649 -10.620  1.00  0.00           S
ATOM     78  CE  MET    10      -6.603   0.179  -9.293  1.00  0.00           C
ATOM     79  N   LEU    11      -3.991  -2.124 -15.047  1.00  0.00
ATOM     80  CA  LEU    11      -3.261  -2.533 -16.235  1.00  0.00
ATOM     81  CB  LEU    11      -1.841  -2.948 -15.844  1.00  0.00
ATOM     82  C   LEU    11      -4.024  -3.658 -16.937  1.00  0.00
ATOM     83  O   LEU    11      -4.213  -3.623 -18.151  1.00  0.00
ATOM     84  CG  LEU    11      -0.981  -1.901 -15.129  1.00  0.00           C
ATOM     85  CD1 LEU    11       0.280  -2.539 -14.566  1.00  0.00           C
ATOM     86  CD2 LEU    11      -0.630  -0.741 -16.052  1.00  0.00           C
ATOM     87  N   ASN    12      -4.443  -4.631 -16.140  1.00  0.00
ATOM     88  CA  ASN    12      -5.183  -5.765 -16.669  1.00  0.00
ATOM     89  CB  ASN    12      -5.506  -6.738 -15.534  1.00  0.00
ATOM     90  C   ASN    12      -6.440  -5.262 -17.382  1.00  0.00
ATOM     91  O   ASN    12      -6.738  -5.685 -18.497  1.00  0.00
ATOM     92  CG  ASN    12      -4.267  -7.283 -14.853  1.00  0.00           C
ATOM     93  OD1 ASN    12      -3.191  -6.693 -14.942  1.00  0.00           O
ATOM     94  ND2 ASN    12      -4.411  -8.410 -14.165  1.00  0.00           N
ATOM     95  N   GLN    13      -7.144  -4.364 -16.707  1.00  0.00
ATOM     96  CA  GLN    13      -8.362  -3.798 -17.261  1.00  0.00
ATOM     97  CB  GLN    13      -8.975  -2.821 -16.256  1.00  0.00
ATOM     98  C   GLN    13      -8.047  -3.133 -18.603  1.00  0.00
ATOM     99  O   GLN    13      -8.755  -3.341 -19.586  1.00  0.00
ATOM    100  CG  GLN    13     -10.314  -2.234 -16.700  1.00  0.00           C
ATOM    101  CD  GLN    13     -10.923  -1.315 -15.659  1.00  0.00           C
ATOM    102  OE1 GLN    13     -10.348  -1.098 -14.593  1.00  0.00           O
ATOM    103  NE2 GLN    13     -12.095  -0.769 -15.964  1.00  0.00           N
ATOM    104  N   PRO    14      -6.981  -2.345 -18.600  1.00  0.00
ATOM    105  CA  PRO    14      -6.562  -1.648 -19.804  1.00  0.00
ATOM    106  CB  PRO    14      -5.330  -0.794 -19.497  1.00  0.00
ATOM    107  C   PRO    14      -6.302  -2.668 -20.915  1.00  0.00
ATOM    108  O   PRO    14      -6.758  -2.492 -22.043  1.00  0.00
ATOM    109  CG  PRO    14      -5.524  -0.554 -18.008  1.00  0.00           C
ATOM    110  CD  PRO    14      -6.192  -1.774 -17.444  1.00  0.00           C
ATOM    111  N   SER    15      -5.570  -3.712 -20.555  1.00  0.00
ATOM    112  CA  SER    15      -5.244  -4.761 -21.507  1.00  0.00
ATOM    113  CB  SER    15      -4.367  -5.814 -20.827  1.00  0.00
ATOM    114  C   SER    15      -6.538  -5.354 -22.070  1.00  0.00
ATOM    115  O   SER    15      -6.671  -5.527 -23.280  1.00  0.00
ATOM    116  OG  SER    15      -4.030  -6.856 -21.738  1.00  0.00           O
ATOM    117  N   ARG    16      -7.458  -5.648 -21.163  1.00  0.00
ATOM    118  CA  ARG    16      -8.737  -6.218 -21.553  1.00  0.00
ATOM    119  CB  ARG    16      -9.580  -6.482 -20.305  1.00  0.00
ATOM    120  C   ARG    16      -9.432  -5.275 -22.538  1.00  0.00
ATOM    121  O   ARG    16      -9.931  -5.711 -23.574  1.00  0.00
ATOM    122  CG  ARG    16      -8.978  -7.483 -19.324  1.00 10.00           C
ATOM    123  CD  ARG    16      -9.087  -8.927 -19.807  1.00 10.00           C
ATOM    124  NE  ARG    16     -10.469  -9.413 -19.812  1.00 10.00           N
ATOM    125  CZ  ARG    16     -11.280  -9.461 -20.871  1.00 10.00           C
ATOM    126  NH1 ARG    16     -10.897  -9.056 -22.082  1.00 10.00           N
ATOM    127  NH2 ARG    16     -12.512  -9.930 -20.716  1.00 10.00           N
ATOM    128  N   THR    17      -9.442  -4.000 -22.179  1.00  0.00
ATOM    129  CA  THR    17     -10.067  -2.991 -23.017  1.00  0.00
ATOM    130  CB  THR    17      -9.955  -1.623 -22.341  1.00  0.00
ATOM    131  C   THR    17      -9.418  -3.011 -24.403  1.00  0.00
ATOM    132  O   THR    17     -10.112  -3.007 -25.418  1.00  0.00
ATOM    133  OG1 THR    17      -8.584  -1.261 -22.154  1.00  0.00           O
ATOM    134  CG2 THR    17     -10.678  -1.596 -21.004  1.00  0.00           C
ATOM    135  N   TRP    18      -8.093  -3.034 -24.400  1.00  0.00
ATOM    136  CA  TRP    18      -7.342  -3.056 -25.643  1.00  0.00
ATOM    137  CB  TRP    18      -5.844  -3.047 -25.335  1.00  0.00
ATOM    138  C   TRP    18      -7.761  -4.276 -26.467  1.00  0.00
ATOM    139  O   TRP    18      -8.021  -4.163 -27.663  1.00  0.00
ATOM    140  CG  TRP    18      -5.388  -1.843 -24.573  1.00  0.00           C
ATOM    141  CD1 TRP    18      -5.324  -1.712 -23.216  1.00  0.00           C
ATOM    142  CD2 TRP    18      -4.938  -0.595 -25.117  1.00  0.00           C
ATOM    143  NE1 TRP    18      -4.862  -0.461 -22.882  1.00  0.00           N
ATOM    144  CE2 TRP    18      -4.618   0.243 -24.032  1.00  0.00           C
ATOM    145  CE3 TRP    18      -4.775  -0.108 -26.418  1.00  0.00           C
ATOM    146  CZ2 TRP    18      -4.145   1.541 -24.205  1.00  0.00           C
ATOM    147  CZ3 TRP    18      -4.305   1.182 -26.589  1.00  0.00           C
ATOM    148  CH2 TRP    18      -3.995   1.992 -25.489  1.00  0.00           C
ATOM    149  N   VAL    19      -7.813  -5.415 -25.792  1.00  0.00
ATOM    150  CA  VAL    19      -8.197  -6.655 -26.445  1.00  0.00
ATOM    151  CB  VAL    19      -8.138  -7.803 -25.436  1.00  0.00
ATOM    152  C   VAL    19      -9.587  -6.493 -27.064  1.00  0.00
ATOM    153  O   VAL    19      -9.803  -6.858 -28.218  1.00  0.00
ATOM    154  CG1 VAL    19      -8.458  -9.140 -26.140  1.00  0.00           C
ATOM    155  CG2 VAL    19      -6.780  -7.889 -24.759  1.00  0.00           C
ATOM    156  N   TYR    20     -10.493  -5.943 -26.268  1.00  0.00
ATOM    157  CA  TYR    20     -11.856  -5.727 -26.723  1.00  0.00
ATOM    158  CB  TYR    20     -12.680  -5.113 -25.590  1.00  0.00
ATOM    159  C   TYR    20     -11.839  -4.848 -27.975  1.00  0.00
ATOM    160  O   TYR    20     -12.504  -5.155 -28.963  1.00  0.00
ATOM    161  CG  TYR    20     -12.783  -5.980 -24.357  1.00  0.00           C
ATOM    162  CD1 TYR    20     -11.878  -5.838 -23.312  1.00  0.00           C
ATOM    163  CD2 TYR    20     -13.778  -6.941 -24.236  1.00  0.00           C
ATOM    164  CE1 TYR    20     -11.961  -6.629 -22.182  1.00  0.00           C
ATOM    165  CE2 TYR    20     -13.869  -7.738 -23.110  1.00  0.00           C
ATOM    166  CZ  TYR    20     -12.958  -7.577 -22.086  1.00  0.00           C
ATOM    167  OH  TYR    20     -13.045  -8.368 -20.963  1.00  0.00           O
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
ATOM      5  CB  ALA A   1     -33.586   2.364  -0.139  1.00 30.00           C
ATOM      6  N   ALA A   2     -31.529  -0.207  -0.287  1.00 30.00           N
ATOM      7  CA  ALA A   2     -30.265  -0.830   0.121  1.00 30.00           C
ATOM      8  C   ALA A   2     -29.092   0.104  -0.268  1.00 30.00           C
ATOM      9  O   ALA A   2     -29.024   0.723  -1.359  1.00 30.00           O
ATOM     10  CB  ALA A   2     -29.837  -2.178  -0.470  1.00 30.00           C
ATOM     11  N   ALA A   3     -28.227   0.361   0.480  1.00 30.00           N
ATOM     12  CA  ALA A   3     -26.960   0.830  -0.121  1.00 30.00           C
ATOM     13  C   ALA A   3     -25.751   0.227   0.567  1.00 30.00           C
ATOM     14  O   ALA A   3     -25.788   0.205   1.784  1.00 30.00           O
ATOM     15  CB  ALA A   3     -26.846   2.364  -0.139  1.00 30.00           C
ATOM     16  N   ALA A   4     -24.789  -0.207  -0.287  1.00 30.00           N
ATOM     17  CA  ALA A   4     -23.525  -0.830   0.121  1.00 30.00           C
ATOM     18  C   ALA A   4     -22.352   0.104  -0.268  1.00 30.00           C
ATOM     19  O   ALA A   4     -22.284   0.723  -1.359  1.00 30.00           O
ATOM     20  CB  ALA A   4     -23.097  -2.178  -0.470  1.00 30.00           C
ATOM     21  N   ALA A   5     -21.487   0.361   0.480  1.00 30.00           N
ATOM     22  CA  ALA A   5     -20.220   0.830  -0.121  1.00 30.00           C
ATOM     23  C   ALA A   5     -19.011   0.227   0.567  1.00 30.00           C
ATOM     24  O   ALA A   5     -19.048   0.205   1.784  1.00 30.00           O
ATOM     25  CB  ALA A   5     -20.106   2.364  -0.139  1.00 30.00           C
ATOM     26  N   ALA A   6     -18.049  -0.207  -0.287  1.00 30.00           N
ATOM     27  CA  ALA A   6     -16.785  -0.830   0.121  1.00 30.00           C
ATOM     28  C   ALA A   6     -15.612   0.104  -0.268  1.00 30.00           C
ATOM     29  O   ALA A   6     -15.544   0.723  -1.359  1.00 30.00           O
ATOM     30  CB  ALA A   6     -16.357  -2.178  -0.470  1.00 30.00           C
ATOM     31  N   ALA A   7     -14.747   0.361   0.480  1.00 30.00           N
ATOM     32  CA  ALA A   7     -13.480   0.830  -0.121  1.00 30.00           C
ATOM     33  C   ALA A   7     -12.271   0.227   0.567  1.00 30.00           C
ATOM     34  O   ALA A   7     -12.308   0.205   1.784  1.00 30.00           O
ATOM     35  CB  ALA A   7     -13.366   2.364  -0.139  1.00 30.00           C
ATOM     36  N   ALA A   8     -11.309  -0.207  -0.287  1.00 30.00           N
ATOM     37  CA  ALA A   8     -10.045  -0.830   0.121  1.00 30.00           C
ATOM     38  C   ALA A   8      -8.872   0.104  -0.268  1.00 30.00           C
ATOM     39  O   ALA A   8      -8.804   0.723  -1.359  1.00 30.00           O
ATOM     40  CB  ALA A   8      -9.617  -2.178  -0.470  1.00 30.00           C
ATOM     41  N   ALA A   9      -8.007   0.361   0.480  1.00 30.00           N
ATOM     42  CA  ALA A   9      -6.740   0.830  -0.121  1.00 30.00           C
ATOM     43  C   ALA A   9      -5.531   0.227   0.567  1.00 30.00           C
ATOM     44  O   ALA A   9      -5.568   0.205   1.784  1.00 30.00           O
ATOM     45  CB  ALA A   9      -6.626   2.364  -0.139  1.00 30.00           C
ATOM     46  N   ALA A  10      -4.569  -0.207  -0.287  1.00 30.00           N
ATOM     47  CA  ALA A  10      -3.305  -0.830   0.121  1.00 30.00           C
ATOM     48  C   ALA A  10      -2.132   0.104  -0.268  1.00 30.00           C
ATOM     49  O   ALA A  10      -2.064   0.723  -1.359  1.00 30.00           O
ATOM     50  CB  ALA A  10      -2.877  -2.178  -0.470  1.00 30.00           C
ATOM     51  N   ALA A  11      -1.267   0.361   0.480  1.00 30.00           N
ATOM     52  CA  ALA A  11       0.000   0.830  -0.121  1.00 30.00           C
ATOM     53  C   ALA A  11       1.209   0.227   0.567  1.00 30.00           C
ATOM     54  O   ALA A  11       1.172   0.205   1.784  1.00 30.00           O
ATOM     55  CB  ALA A  11       0.114   2.364  -0.139  1.00 30.00           C
ATOM     56  N   ALA A  12       2.171  -0.207  -0.287  1.00 30.00           N
ATOM     57  CA  ALA A  12       3.435  -0.830   0.121  1.00 30.00           C
ATOM     58  C   ALA A  12       4.608   0.104  -0.268  1.00 30.00           C
ATOM     59  O   ALA A  12       4.676   0.723  -1.359  1.00 30.00           O
ATOM     60  CB  ALA A  12       3.863  -2.178  -0.470  1.00 30.00           C
ATOM     61  N   ALA A  13       5.473   0.361   0.480  1.00 30.00           N
ATOM     62  CA  ALA A  13       6.740   0.830  -0.121  1.00 30.00           C
ATOM     63  C   ALA A  13       7.949   0.227   0.567  1.00 30.00           C
ATOM     64  O   ALA A  13       7.912   0.205   1.784  1.00 30.00           O
ATOM     65  CB  ALA A  13       6.854   2.364  -0.139  1.00 30.00           C
ATOM     66  N   ALA A  14       8.911  -0.207  -0.287  1.00 30.00           N
ATOM     67  CA  ALA A  14      10.175  -0.830   0.121  1.00 30.00           C
ATOM     68  C   ALA A  14      11.348   0.104  -0.268  1.00 30.00           C
ATOM     69  O   ALA A  14      11.416   0.723  -1.359  1.00 30.00           O
ATOM     70  CB  ALA A  14      10.603  -2.178  -0.470  1.00 30.00           C
ATOM     71  N   ALA A  15      12.213   0.361   0.480  1.00 30.00           N
ATOM     72  CA  ALA A  15      13.480   0.830  -0.121  1.00 30.00           C
ATOM     73  C   ALA A  15      14.689   0.227   0.567  1.00 30.00           C
ATOM     74  O   ALA A  15      14.652   0.205   1.784  1.00 30.00           O
ATOM     75  CB  ALA A  15      13.594   2.364  -0.139  1.00 30.00           C
ATOM     76  N   ALA A  16      15.651  -0.207  -0.287  1.00 30.00           N
ATOM     77  CA  ALA A  16      16.915  -0.830   0.121  1.00 30.00           C
ATOM     78  C   ALA A  16      18.088   0.104  -0.268  1.00 30.00           C
ATOM     79  O   ALA A  16      18.156   0.723  -1.359  1.00 30.00           O
ATOM     80  CB  ALA A  16      17.343  -2.178  -0.470  1.00 30.00           C
ATOM     81  N   ALA A  17      18.953   0.361   0.480  1.00 30.00           N
ATOM     82  CA  ALA A  17      20.220   0.830  -0.121  1.00 30.00           C
ATOM     83  C   ALA A  17      21.429   0.227   0.567  1.00 30.00           C
ATOM     84  O   ALA A  17      21.392   0.205   1.784  1.00 30.00           O
ATOM     85  CB  ALA A  17      20.334   2.364  -0.139  1.00 30.00           C
ATOM     86  N   ALA A  18      22.391  -0.207  -0.287  1.00 30.00           N
ATOM     87  CA  ALA A  18      23.655  -0.830   0.121  1.00 30.00           C
ATOM     88  C   ALA A  18      24.828   0.104  -0.268  1.00 30.00           C
ATOM     89  O   ALA A  18      24.896   0.723  -1.359  1.00 30.00           O
ATOM     90  CB  ALA A  18      24.083  -2.178  -0.470  1.00 30.00           C
ATOM     91  N   ALA A  19      25.693   0.361   0.480  1.00 30.00           N
ATOM     92  CA  ALA A  19      26.960   0.830  -0.121  1.00 30.00           C
ATOM     93  C   ALA A  19      28.169   0.227   0.567  1.00 30.00           C
ATOM     94  O   ALA A  19      28.132   0.205   1.784  1.00 30.00           O
ATOM     95  CB  ALA A  19      27.074   2.364  -0.139  1.00 30.00           C
ATOM     96  N   ALA A  20      29.131  -0.207  -0.287  1.00 30.00           N
ATOM     97  CA  ALA A  20      30.395  -0.830   0.121  1.00 30.00           C
ATOM     98  C   ALA A  20      31.568   0.104  -0.268  1.00 30.00           C
ATOM     99  O   ALA A  20      31.636   0.723  -1.359  1.00 30.00           O
ATOM    100  CB  ALA A  20      30.823  -2.178  -0.470  1.00 30.00           C
ATOM    101  N   ALA A  21      32.433   0.361   0.480  1.00 30.00           N
ATOM    102  CA  ALA A  21      33.700   0.830  -0.121  1.00 30.00           C
ATOM    103  C   ALA A  21      34.909   0.227   0.567  1.00 30.00           C
ATOM    104  O   ALA A  21      34.872   0.205   1.784  1.00 30.00           O
ATOM    105  CB  ALA A  21      33.814   2.364  -0.139  1.00 30.00           C
ATOM    106  N   ALA A  22      35.871  -0.207  -0.287  1.00 30.00           N
ATOM    107  CA  ALA A  22      37.135  -0.830   0.121  1.00 30.00           C
ATOM    108  C   ALA A  22      38.308   0.104  -0.268  1.00 30.00           C
ATOM    109  O   ALA A  22      38.376   0.723  -1.359  1.00 30.00           O
ATOM    110  CB  ALA A  22      37.563  -2.178  -0.470  1.00 30.00           C
TER
END
"""

def exercise_00(prefix="exercise_00"):
  "Build poly-ALA alpha helix"
  ph = ssb.make_ss_structure_from_seq(
    alpha_helix_template, "".join(["A"]*12))
  ph.write_pdb_file(file_name="%s_result.pdb"%prefix)
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
  assert d2.min_max_mean().as_tuple()[2]<0.05

def exercise_01(prefix="exercise_01"):
  "Build poly-ALA beta strand"
  ph = ssb.make_ss_structure_from_seq(
    beta_strand_template, "".join(["A"]*22))
  ph.write_pdb_file(file_name="%s_result.pdb"%prefix)
  of=open("%s_answer.pdb"%prefix, "w")
  print >> of, beta_strand_answer
  of.close()
  sites_1 = ph.atoms().extract_xyz()
  sites_2 = iotbx.pdb.input(source_info=None,
    lines=beta_strand_answer).construct_hierarchy().atoms().extract_xyz()
  assert sites_1.size()==sites_2.size()

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
  except Sorry: pass
  else: raise Exception_expected

def exercise_ss_structure_from_seq():
  pdb_inp = iotbx.pdb.input(source_info=None, lines=correct_answer)
  correct_h = pdb_inp.construct_hierarchy()
  test_h = ssb.make_ss_structure_from_seq(alpha_pdb_str, "ACEDGFIHKMLNQPSRTWVY")

  # why this assert fails:
  #assert correct_h.is_similar_hierarchy(other=test_h)

  # but this is not:
  #assert test_h.as_str() == correct_h.as_str()

  assert approx_equal(test_h.atoms().extract_xyz(),
                      correct_h.atoms().extract_xyz(), eps=0.002)
  try: ssb.make_ss_structure_from_seq(alpha_pdb_str, "")
  except Sorry: pass
  else: raise Exception_expected

def exercise():
  exercise_00()
  exercise_01()
  exercise_r_t_matrices()
  exercise_ss_structure_from_seq()

if (__name__ == "__main__"):
  exercise()
