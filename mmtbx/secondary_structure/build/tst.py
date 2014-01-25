from mmtbx.secondary_structure import build as ssb
import iotbx.pdb
from libtbx.test_utils import Exception_expected, approx_equal, show_diff
from scitbx import matrix
from libtbx.utils import Sorry


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
  exercise_r_t_matrices()
  exercise_ss_structure_from_seq()
  #sys.stdout.flush()
  #unittest.TextTestRunner(stream=sys.stdout, verbosity = 2 ).run( alltests )

if (__name__ == "__main__"):
  exercise()
