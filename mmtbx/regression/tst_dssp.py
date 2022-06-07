
from __future__ import absolute_import, division, print_function
from libtbx.test_utils import show_diff
from libtbx import easy_run
import libtbx.load_env
import os

def exercise_basic():
  # extracted from 1ywf
  pdb_in = """\
ATOM      1  N   ASP A 181      27.197  36.852   5.496  1.00 33.48           N
ATOM      2  CA  ASP A 181      27.251  36.418   4.117  1.00 34.90           C
ATOM      3  C   ASP A 181      26.472  37.390   3.236  1.00 34.37           C
ATOM      4  O   ASP A 181      25.509  38.039   3.694  1.00 33.69           O
ATOM      5  CB  ASP A 181      26.657  34.998   3.964  1.00 36.77           C
ATOM      6  N   ARG A 182      26.880  37.461   1.974  1.00 33.47           N
ATOM      7  CA  ARG A 182      26.256  38.318   1.002  1.00 33.32           C
ATOM      8  C   ARG A 182      24.744  38.107   0.932  1.00 31.55           C
ATOM      9  O   ARG A 182      23.984  39.056   0.794  1.00 29.58           O
ATOM     10  CB  ARG A 182      26.877  38.089  -0.369  1.00 34.00           C
ATOM     11  N   ASP A 183      24.289  36.868   1.010  1.00 31.07           N
ATOM     12  CA  ASP A 183      22.878  36.632   0.774  1.00 31.07           C
ATOM     13  C   ASP A 183      22.032  37.273   1.878  1.00 29.65           C
ATOM     14  O   ASP A 183      20.955  37.773   1.618  1.00 27.98           O
ATOM     15  CB  ASP A 183      22.592  35.148   0.636  1.00 32.71           C
ATOM     16  N   VAL A 184      22.537  37.250   3.102  1.00 28.40           N
ATOM     17  CA  VAL A 184      21.844  37.855   4.221  1.00 28.15           C
ATOM     18  C   VAL A 184      21.822  39.403   4.109  1.00 26.40           C
ATOM     19  O   VAL A 184      20.785  40.046   4.308  1.00 25.30           O
ATOM     20  CB  VAL A 184      22.533  37.460   5.501  1.00 28.88           C
ATOM     21  N   ILE A 185      22.961  39.970   3.749  1.00 25.19           N
ATOM     22  CA  ILE A 185      23.099  41.414   3.509  1.00 23.33           C
ATOM     23  C   ILE A 185      22.123  41.852   2.412  1.00 23.71           C
ATOM     24  O   ILE A 185      21.376  42.811   2.561  1.00 22.57           O
ATOM     25  CB  ILE A 185      24.556  41.728   3.091  1.00 23.73           C
ATOM     26  N   VAL A 186      22.112  41.120   1.311  1.00 23.54           N
ATOM     27  CA  VAL A 186      21.221  41.455   0.194  1.00 24.20           C
ATOM     28  C   VAL A 186      19.738  41.300   0.563  1.00 23.20           C
ATOM     29  O   VAL A 186      18.906  42.075   0.105  1.00 24.19           O
ATOM     30  CB  VAL A 186      21.584  40.632  -1.092  1.00 24.98           C
ATOM     31  N   ALA A 187      19.403  40.316   1.385  1.00 23.05           N
ATOM     32  CA  ALA A 187      18.035  40.114   1.854  1.00 23.26           C
ATOM     33  C   ALA A 187      17.523  41.382   2.547  1.00 22.44           C
ATOM     34  O   ALA A 187      16.427  41.902   2.256  1.00 21.80           O
ATOM     35  CB  ALA A 187      17.988  38.991   2.765  1.00 23.85           C
ATOM     36  N   ASP A 188      18.365  41.932   3.397  1.00 21.03           N
ATOM     37  CA  ASP A 188      17.956  43.118   4.149  1.00 21.44           C
ATOM     38  C   ASP A 188      17.939  44.349   3.264  1.00 20.78           C
ATOM     39  O   ASP A 188      17.042  45.182   3.399  1.00 22.01           O
ATOM     40  CB  ASP A 188      18.867  43.342   5.384  1.00 21.37           C
ATOM     41  N   TYR A 189      18.896  44.462   2.354  1.00 20.83           N
ATOM     42  CA  TYR A 189      18.982  45.542   1.424  1.00 20.65           C
ATOM     43  C   TYR A 189      17.722  45.582   0.535  1.00 21.46           C
ATOM     44  O   TYR A 189      17.134  46.647   0.250  1.00 20.62           O
ATOM     45  CB  TYR A 189      20.232  45.356   0.573  1.00 21.85           C
ATOM     46  N   LEU A 190      17.306  44.409   0.093  1.00 22.68           N
ATOM     47  CA  LEU A 190      16.176  44.366  -0.846  1.00 22.91           C
ATOM     48  C   LEU A 190      14.848  44.614  -0.155  1.00 22.51           C
ATOM     49  O   LEU A 190      13.859  44.905  -0.830  1.00 22.82           O
ATOM     50  CB  LEU A 190      16.116  43.032  -1.594  1.00 22.51           C
ATOM     51  N   ARG A 191      14.793  44.581   1.175  1.00 22.46           N
ATOM     52  CA  ARG A 191      13.551  44.901   1.906  1.00 23.66           C
ATOM     53  C   ARG A 191      13.025  46.291   1.588  1.00 22.52           C
ATOM     54  O   ARG A 191      11.824  46.554   1.794  1.00 21.76           O
ATOM     55  CB  ARG A 191      13.735  44.792   3.436  1.00 24.11           C
ATOM     56  N   SER A 192      13.878  47.175   1.057  1.00 22.09           N
ATOM     57  CA  SER A 192      13.454  48.520   0.682  1.00 21.85           C
ATOM     58  C   SER A 192      12.449  48.512  -0.485  1.00 22.55           C
ATOM     59  O   SER A 192      11.624  49.409  -0.607  1.00 22.52           O
ATOM     60  CB  SER A 192      14.639  49.422   0.322  1.00 22.37           C
ATOM     61  N   ASN A 193      12.467  47.441  -1.269  1.00 22.93           N
ATOM     62  CA  ASN A 193      11.435  47.276  -2.282  1.00 23.18           C
ATOM     63  C   ASN A 193      10.046  47.185  -1.736  1.00 22.19           C
ATOM     64  O   ASN A 193       9.094  47.629  -2.388  1.00 23.39           O
ATOM     65  CB  ASN A 193      11.770  46.077  -3.193  1.00 24.30           C
ATOM     66  N   ASP A 194       9.878  46.683  -0.515  1.00 23.18           N
ATOM     67  CA  ASP A 194       8.554  46.586   0.069  1.00 23.55           C
ATOM     68  C   ASP A 194       7.939  47.961   0.370  1.00 24.28           C
ATOM     69  O   ASP A 194       6.720  48.105   0.555  1.00 23.59           O
ATOM     70  CB  ASP A 194       8.605  45.788   1.353  1.00 22.75           C
ATOM     71  N   SER A 195       8.809  48.958   0.499  1.00 23.92           N
ATOM     72  CA  SER A 195       8.409  50.341   0.846  1.00 24.19           C
ATOM     73  C   SER A 195       8.150  51.237  -0.368  1.00 23.29           C
ATOM     74  O   SER A 195       7.884  52.429  -0.240  1.00 23.19           O
ATOM     75  CB  SER A 195       9.521  50.937   1.693  1.00 24.32           C
ATOM     76  N   VAL A 196       8.253  50.699  -1.568  1.00 23.63           N
ATOM     77  CA  VAL A 196       8.152  51.513  -2.768  1.00 23.74           C
ATOM     78  C   VAL A 196       6.773  52.193  -2.856  1.00 23.84           C
ATOM     79  O   VAL A 196       6.718  53.397  -3.105  1.00 23.11           O
ATOM     80  CB  VAL A 196       8.465  50.710  -4.027  1.00 24.47           C
ATOM     81  N   PRO A 197       5.665  51.486  -2.598  1.00 24.39           N
ATOM     82  CA  PRO A 197       4.365  52.177  -2.602  1.00 26.12           C
ATOM     83  C   PRO A 197       4.278  53.337  -1.596  1.00 26.03           C
ATOM     84  O   PRO A 197       3.742  54.399  -1.930  1.00 25.76           O
ATOM     85  CB  PRO A 197       3.363  51.049  -2.246  1.00 26.62           C
ATOM     86  N   GLN A 198       4.859  53.168  -0.408  1.00 25.67           N
ATOM     87  CA  GLN A 198       4.860  54.206   0.603  1.00 26.81           C
ATOM     88  C   GLN A 198       5.732  55.384   0.158  1.00 25.59           C
ATOM     89  O   GLN A 198       5.341  56.545   0.339  1.00 25.72           O
ATOM     90  CB  GLN A 198       5.290  53.638   1.967  1.00 27.38           C
ATOM     91  N   LEU A 199       6.903  55.119  -0.416  1.00 25.48           N
ATOM     92  CA  LEU A 199       7.726  56.192  -0.941  1.00 25.93           C
ATOM     93  C   LEU A 199       6.996  56.972  -2.047  1.00 26.39           C
ATOM     94  O   LEU A 199       7.020  58.180  -2.064  1.00 25.38           O
ATOM     95  CB  LEU A 199       9.033  55.633  -1.490  1.00 25.66           C
ATOM     96  N   ARG A 200       6.361  56.258  -2.980  1.00 27.31           N
ATOM     97  CA  ARG A 200       5.576  56.913  -3.993  1.00 28.53           C
ATOM     98  C   ARG A 200       4.520  57.823  -3.397  1.00 27.54           C
ATOM     99  O   ARG A 200       4.397  58.949  -3.851  1.00 27.50           O
ATOM    100  CB  ARG A 200       4.933  55.879  -4.899  1.00 30.38           C
ATOM    101  N   ALA A 201       3.790  57.365  -2.357  1.00 26.90           N
ATOM    102  CA  ALA A 201       2.764  58.200  -1.713  1.00 26.49           C
ATOM    103  C   ALA A 201       3.406  59.407  -1.045  1.00 26.59           C
ATOM    104  O   ALA A 201       2.866  60.516  -1.082  1.00 26.58           O
ATOM    105  CB  ALA A 201       1.959  57.412  -0.715  1.00 25.11           C
ATOM    106  N   ARG A 202       4.566  59.205  -0.419  1.00 25.66           N
ATOM    107  CA  ARG A 202       5.245  60.296   0.240  1.00 26.93           C
ATOM    108  C   ARG A 202       5.676  61.346  -0.767  1.00 26.93           C
ATOM    109  O   ARG A 202       5.555  62.541  -0.489  1.00 25.79           O
ATOM    110  CB  ARG A 202       6.493  59.779   0.996  1.00 28.25           C
ATOM    111  N   ILE A 203       6.154  60.912  -1.931  1.00 26.99           N
ATOM    112  CA  ILE A 203       6.611  61.848  -2.965  1.00 27.49           C
ATOM    113  C   ILE A 203       5.430  62.674  -3.480  1.00 28.29           C
ATOM    114  O   ILE A 203       5.548  63.905  -3.624  1.00 27.82           O
ATOM    115  CB  ILE A 203       7.322  61.125  -4.075  1.00 28.09           C
ATOM    116  N   SER A 204       4.288  62.025  -3.678  1.00 27.96           N
ATOM    117  CA  SER A 204       3.119  62.736  -4.184  1.00 28.04           C
ATOM    118  C   SER A 204       2.683  63.793  -3.199  1.00 28.16           C
ATOM    119  O   SER A 204       2.311  64.910  -3.605  1.00 28.25           O
ATOM    120  CB  SER A 204       1.962  61.780  -4.504  1.00 27.64           C
ATOM    121  N   GLU A 205       2.756  63.477  -1.904  1.00 27.27           N
ATOM    122  CA  GLU A 205       2.419  64.414  -0.849  1.00 28.81           C
ATOM    123  C   GLU A 205       3.430  65.594  -0.782  1.00 28.51           C
ATOM    124  O   GLU A 205       3.035  66.755  -0.615  1.00 29.81           O
ATOM    125  CB  GLU A 205       2.370  63.667   0.461  1.00 29.74           C
ATOM    126  N   MET A 206       4.709  65.295  -0.984  1.00 28.00           N
ATOM    127  CA  MET A 206       5.726  66.347  -1.029  1.00 27.28           C
ATOM    128  C   MET A 206       5.494  67.315  -2.173  1.00 27.94           C
ATOM    129  O   MET A 206       5.672  68.524  -2.010  1.00 27.36           O
ATOM    130  CB  MET A 206       7.152  65.714  -1.123  1.00 27.25           C
ATOM    131  N   ILE A 207       5.092  66.772  -3.319  1.00 29.22           N
ATOM    132  CA  ILE A 207       4.824  67.562  -4.510  1.00 29.58           C
ATOM    133  C   ILE A 207       3.606  68.471  -4.277  1.00 30.18           C
ATOM    134  O   ILE A 207       3.601  69.645  -4.680  1.00 28.59           O
ATOM    135  CB  ILE A 207       4.648  66.646  -5.716  1.00 30.41           C
ATOM    136  N   GLN A 208       2.586  67.940  -3.614  1.00 29.64           N
ATOM    137  CA  GLN A 208       1.415  68.721  -3.231  1.00 30.71           C
ATOM    138  C   GLN A 208       1.799  69.877  -2.302  1.00 29.77           C
ATOM    139  O   GLN A 208       1.213  70.964  -2.365  1.00 30.52           O
ATOM    140  CB  GLN A 208       0.345  67.845  -2.558  1.00 30.05           C
ATOM    141  N   GLN A 209       2.776  69.641  -1.418  1.00 29.19           N
ATOM    142  CA  GLN A 209       3.254  70.616  -0.496  1.00 28.05           C
ATOM    143  C   GLN A 209       4.454  71.438  -1.005  1.00 26.36           C
ATOM    144  O   GLN A 209       5.077  72.111  -0.213  1.00 26.13           O
ATOM    145  CB  GLN A 209       3.663  69.913   0.841  1.00 28.70           C
ATOM    146  N   ARG A 210       4.719  71.435  -2.303  1.00 26.63           N
ATOM    147  CA  ARG A 210       5.975  72.001  -2.823  1.00 27.68           C
ATOM    148  C   ARG A 210       6.081  73.522  -2.641  1.00 28.32           C
ATOM    149  O   ARG A 210       7.167  74.068  -2.709  1.00 28.23           O
ATOM    150  CB  ARG A 210       6.178  71.657  -4.275  1.00 27.92           C
END"""
  with open("tst_dssp_1ywf_helix.pdb", "w") as f:
    f.write(pdb_in)
  result = easy_run.fully_buffered("mmtbx.dssp tst_dssp_1ywf_helix.pdb"
    ).raise_if_errors()
  assert ("\n".join(result.stdout_lines) == """\
HELIX    1   1 ASP A  181  ARG A  191  1                                  11
HELIX    3   3 SER A  195  GLN A  209  1                                  15""")

# these exercised depend on files in other SVN trees
def exercise_advanced():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  if (pdb_file is None):
    print("phenix_regression not available, skipping test")
    return
  result = easy_run.fully_buffered("mmtbx.dssp \"%s\"" % pdb_file
    ).raise_if_errors()
  assert not show_diff("\n".join(result.stdout_lines), """\
HELIX    2   2 ASP A   14  THR A   18  5                                   5
HELIX    6   6 ASP A   37  GLY A   48  1                                  12
HELIX    8   8 SER A   57  GLY A   65  1                                   9
HELIX   10  10 ASN A  119  GLN A  132  1                                  14
HELIX   16  16 GLY A  138  ALA A  152  1                                  15
HELIX   19  19 ASP A  165  VAL A  178  1                                  14
HELIX   22  22 ASP A  181  ARG A  191  1                                  11
HELIX   24  24 SER A  195  GLN A  209  1                                  15
HELIX   29  29 ALA A  216  ALA A  225  1                                  10
HELIX   33  33 SER A  228  GLY A  233  1                                   6
HELIX   35  35 ARG A  235  TYR A  250  1                                  16
HELIX   38  38 SER A  252  ALA A  260  1                                   9
HELIX   41  41 SER A  263  LEU A  275  1                                  13
SHEET    1   1 5 ARG A  13  ASP A  14  0
SHEET    2   1 5 LEU A  27  SER A  30 -1  O  ARG A  29   N  ARG A  13
SHEET    3   1 5 VAL A 156  HIS A 159  1  O  VAL A 156   N  PHE A  28
SHEET    4   1 5 ASP A  51  ASP A  54  1  O  ASP A  51   N  LEU A 157
SHEET    5   1 5 ASP A  74  LEU A  77  1  O  ASP A  74   N  VAL A  52""")
  cctbx_p_dir = libtbx.env.find_in_repositories(
    relative_path="cctbx_project",
    test=os.path.isdir)
  if (cctbx_p_dir is not None):
    # mostly beta
    pdb_file_1 = os.path.join(cctbx_p_dir, "mmtbx", "regression", "pdbs", "p9.pdb")
    result = easy_run.fully_buffered("mmtbx.dssp \"%s\"" % pdb_file_1
      ).raise_if_errors()
    # XXX interestingly, this appears to be more accurate than the output of
    # ksdssp, which finds an alpha helix from 111 to 116.  both the PDB and
    # ksdssp have 131-137 as a continuous strand, but residue 135 actually
    # breaks the hydrogen bonding pattern, and the hydrogen bonds cannot be
    # accurately recovered from those annotations.
    assert not show_diff("\n".join(result.stdout_lines), """\
HELIX   10  10 PRO    105  VAL    109  5                                   5
HELIX   13  13 ALA    113  LEU    117  5                                   5
SHEET    1   1 2 TYR    11  GLU    13  0
SHEET    2   1 2 GLN    70  GLU    72 -1  O  VAL    71   N  VAL    12
SHEET    1   2 4 TYR    22  ILE    25  0
SHEET    2   2 4 GLU    28  VAL    32 -1  O  CYS    30   N  VAL    23
SHEET    3   2 4 LYS    47  GLY    54 -1  O  VAL    53   N  ARG    31
SHEET    4   2 4 LYS    60  PRO    66 -1  O  LEU    65   N  ALA    48
SHEET    1   3 2 GLU    34  SER    38  0
SHEET    2   3 2 LYS    47  VAL    51 -1  O  VAL    51   N  GLU    34
SHEET    1   4 5 SER    85  VAL    86  0
SHEET    2   4 5 ILE    91  ASP    95 -1  O  GLN    92   N  SER    85
SHEET    3   4 5 GLU    77  ILE    83 -1  O  GLN    82   N  MET    94
SHEET    4   4 5 GLU   122  ILE   128 -1  O  GLN   127   N  GLU    77
SHEET    5   4 5 ARG   131  ILE   134 -1  O  LYS   133   N  TRP   126
SHEET    1   5 2 VAL    90  MET    94  0
SHEET    2   5 2 THR   101  PRO   105 -1  O  VAL   104   N  ILE    91
SHEET    1   6 2 VAL   123  GLU   124  0
SHEET    2   6 2 ARG   136  VAL   137 -1  O  ARG   136   N  GLU   124""")
    # beta barrel
  examples_dir = libtbx.env.find_in_repositories(
    relative_path="phenix_examples",
    test=os.path.isdir)
  if examples_dir is not None:
    pdb_file_2 = os.path.join(examples_dir, "porin-twin", "porin.pdb")
    result = easy_run.fully_buffered("mmtbx.dssp \"%s\"" % pdb_file_2
      ).raise_if_errors()
    assert not show_diff("\n".join(result.stdout_lines), """\
HELIX    5   5 ASP     59  GLY     63  5                                   5
HELIX    7   7 THR     87  VAL     92  1                                   6
HELIX   18  18 ASP    156  VAL    160  5                                   5
HELIX   19  19 ASP    184  ILE    188  5                                   5
SHEET    1   117 ILE     2  TYR    14  0
SHEET    2   117 THR    25  GLU    40 -1  O  VAL    36   N  SER     3
SHEET    3   117 THR    46  ASP    56 -1  O  TRP    55   N  LEU    31
SHEET    4   117 GLN    70  TYR    75 -1  O  SER    74   N  THR    46
SHEET    5   117 VAL    78  GLY    83 -1  O  VAL    82   N  PHE    71
SHEET    6   117 ASN   131  ILE   139 -1  O  THR   136   N  THR    79
SHEET    7   117 VAL   142  ASP   150 -1  O  ASP   150   N  ASN   131
SHEET    8   117 GLU   163  SER   171 -1  O  ASP   169   N  ASN   143
SHEET    9   117 ILE   175  THR   183 -1  O  THR   183   N  PHE   164
SHEET   10   117 ILE   193  LYS   201 -1  O  ALA   199   N  SER   176
SHEET   11   117 GLY   206  ASP   214 -1  O  ASP   214   N  ALA   194
SHEET   12   117 GLN   223  PHE   232 -1  O  ASN   229   N  THR   207
SHEET   13   117 THR   235  ILE   244 -1  O  ASP   243   N  VAL   224
SHEET   14   117 ALA   252  GLN   260 -1  O  ASP   258   N  THR   236
SHEET   15   117 VAL   265  SER   273 -1  O  SER   273   N  TYR   253
SHEET   16   117 THR   279  ASP   288 -1  O  ARG   286   N  LYS   266
SHEET   17   117 TYR     7  VAL    15 -1  O  TYR    14   N  ALA   281""")
  else :
    print("WARNING: phenix_examples not available, some tests skipped")
  # TODO add 1lfh test
  print("OK")

if (__name__ == "__main__"):
  exercise_basic()
  exercise_advanced()
