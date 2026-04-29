"""Test atom selection string usage"""
from __future__ import absolute_import, division, print_function
from iotbx.pdb.atom_selection import selection_string_from_selection
from iotbx.pdb.atom_selection import get_clean_selection_string
from mmtbx.ncs.ncs_search import get_chains_info
from scitbx.array_family import flex
import iotbx.pdb
from libtbx.test_utils import approx_equal
from six.moves import range
from six.moves import zip

test_pdb_1 = '''\
CRYST1  577.812  448.715  468.790  90.00  90.00  90.00 P 1
SCALE1      0.001731  0.000000  0.000000        0.00000
SCALE2      0.000000  0.002229  0.000000        0.00000
SCALE3      0.000000  0.000000  0.002133        0.00000
ATOM      1  CA  LYS A 151      10.766   9.333  12.905  1.00 44.22           C
ATOM      2  CA  LYS A 152      10.117   9.159  11.610  1.00 49.42           C
ATOM      3  CA  LYS A 153       9.099   8.000  11.562  1.00 46.15           C
ATOM      4  CA  LYS A 154       8.000   8.202  11.065  1.00 52.97           C
ATOM      5  CA  LYS A 155      11.146   9.065  10.474  1.00 41.68           C
ATOM      6  CA  LYS A 156      10.547   9.007   9.084  1.00 55.55           C
ATOM      7  CA  LYS A 157      11.545   9.413   8.000  1.00 72.27           C
ATOM      8  CA  LYS A 158      12.277  10.718   8.343  1.00 75.78           C
ATOM      9  CA  LYS A 159      11.349  11.791   8.809  1.00 75.88           C
TER
ATOM    222  CA  LEU X  40      94.618  -5.253  91.582  1.00 87.10           C
ATOM    223  CA  ARG X  41      62.395  51.344  80.786  1.00107.25           C
ATOM    224  CA  ARG X  42      62.395  41.344  80.786  1.00107.25           C
TER
ATOM      1  CA  THR D   1       8.111  11.080  10.645  1.00 20.00           C
ATOM      2  CA  THR D   2       8.000   9.722  10.125  1.00 20.00           C
ATOM      3  CA  THR D   3       8.075   8.694  11.249  1.00 20.00           C
ATOM      4  CA  THR D   4       8.890   8.818  12.163  1.00 20.00           C
ATOM      5  CA  THR D   5       9.101   9.421   9.092  1.00 20.00           C
ATOM      6  CA  THR D   6       9.001  10.343   8.000  1.00 20.00           C
ATOM      7  CA  THR D   7       8.964   8.000   8.565  1.00 20.00           C
END
'''

test_pdb_2 = '''\
CRYST1   41.870   78.240   85.900  90.00  90.00  90.00 P 2 21 21     8
SCALE1      0.023883  0.000000  0.000000        0.00000
SCALE2      0.000000  0.012781  0.000000        0.00000
SCALE3      0.000000  0.000000  0.011641        0.00000
ATOM      1  N   LYS A 151      16.915  16.113 -32.818  1.00 44.22           N
ATOM      2  CA  LYS A 151      16.266  15.939 -34.113  1.00 49.42           C
ATOM      3  C   LYS A 151      15.248  14.780 -34.161  1.00 46.15           C
ATOM      4  O   LYS A 151      14.149  14.982 -34.658  1.00 52.97           O
ATOM      5  CB  LYS A 151      17.295  15.845 -35.249  1.00 41.68           C
ATOM      6  CG  LYS A 151      16.696  15.787 -36.639  1.00 55.55           C
ATOM      7  CD  LYS A 151      17.694  16.193 -37.723  1.00 72.27           C
ATOM      8  CE  LYS A 151      18.426  17.498 -37.380  1.00 75.78           C
ATOM      9  NZ  LYS A 151      17.498  18.571 -36.914  1.00 75.88           N
ATOM     10  N   ARG A 152      15.575  13.584 -33.658  1.00 41.11           N
ATOM     11  CA  ARG A 152      14.545  12.524 -33.623  1.00 39.90           C
ATOM     12  C   ARG A 152      14.775  11.251 -32.810  1.00 33.61           C
ATOM     13  O   ARG A 152      15.911  10.847 -32.533  1.00 34.25           O
ATOM     14  CB  ARG A 152      14.096  12.135 -35.040  1.00 41.02           C
ATOM     15  CG  ARG A 152      15.094  11.299 -35.813  1.00 45.33           C
ATOM     16  CD  ARG A 152      14.593  11.085 -37.231  1.00 44.56           C
ATOM     17  NE  ARG A 152      13.434  10.190 -37.278  1.00 44.10           N
ATOM     18  CZ  ARG A 152      12.734   9.937 -38.378  1.00 49.14           C
ATOM     19  NH1 ARG A 152      13.064  10.532 -39.521  1.00 44.79           N
ATOM     20  NH2 ARG A 152      11.700   9.097 -38.335  1.00 45.30           N
ATOM     21  N   ALA B 153      13.656  10.621 -32.456  1.00 30.81           N
ATOM     22  CA  ALA B 153      13.619   9.397 -31.651  1.00 29.91           C
ATOM     23  C   ALA B 153      14.208   8.219 -32.409  1.00 31.80           C
ATOM     24  O   ALA B 153      14.248   8.241 -33.644  1.00 31.66           O
ATOM     25  CB  ALA B 153      12.176   9.082 -31.250  1.00 29.44           C
ATOM     26  N   PRO B 154      14.636   7.168 -31.680  1.00 30.22           N
ATOM     27  CA  PRO B 154      15.163   5.989 -32.375  1.00 28.83           C
ATOM     28  C   PRO B 154      14.085   5.348 -33.254  1.00 33.21           C
ATOM     29  O   PRO B 154      12.885   5.472 -32.970  1.00 26.91           O
ATOM     30  CB  PRO B 154      15.534   5.028 -31.234  1.00 27.52           C
ATOM     31  CG  PRO B 154      15.512   5.870 -29.971  1.00 30.59           C
ATOM     32  CD  PRO B 154      14.507   6.943 -30.231  1.00 27.77           C
ATOM     33  N   TYR B 155      14.519   4.680 -34.318  1.00 28.01           N
ATOM     34  CA  TYR B 155      13.603   3.970 -35.194  1.00 26.28           C
ATOM     35  C   TYR B 155      14.387   2.846 -35.845  1.00 32.08           C
ATOM     36  O   TYR B 155      15.623   2.936 -35.993  1.00 29.48           O
ATOM     37  CB  TYR B 155      13.028   4.912 -36.264  1.00 29.02           C
ATOM     38  CG  TYR B 155      14.088   5.528 -37.164  1.00 36.38           C
ATOM     39  CD1 TYR B 155      14.789   6.670 -36.769  1.00 29.69           C
ATOM     40  CD2 TYR B 155      14.403   4.957 -38.401  1.00 38.19           C
ATOM     41  CE1 TYR B 155      15.760   7.237 -37.584  1.00 37.04           C
ATOM     42  CE2 TYR B 155      15.376   5.518 -39.226  1.00 43.21           C
ATOM     43  CZ  TYR B 155      16.051   6.654 -38.804  1.00 41.97           C
END
'''

test_pdb_3 = '''\
CRYST1  203.106   83.279  178.234  90.00 106.67  90.00 C 1 2 1      12
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.004924  0.000000  0.001474        0.00000
SCALE2      0.000000  0.012008  0.000000        0.00000
SCALE3      0.000000  0.000000  0.005857        0.00000
ATOM    313  N   CYS H  47      85.603 -27.032   6.791  1.00 42.51           N
ATOM    314  CA  CYS H  47      84.850 -28.275   6.785  1.00 44.04           C
ATOM    315  C   CYS H  47      83.442 -28.089   6.261  1.00 42.34           C
ATOM    316  O   CYS H  47      83.056 -26.989   5.881  1.00 42.28           O
ATOM    317  CB  CYS H  47      84.827 -28.861   8.204  1.00 47.91           C
ATOM    318  SG  CYS H  47      86.496 -29.154   8.879  1.00 51.94           S
ATOM    319  N   ARG H  48      82.680 -29.174   6.224  1.00 40.63           N
ATOM    320  CA  ARG H  48      81.318 -29.102   5.747  1.00 40.99           C
ATOM    321  C   ARG H  48      80.335 -28.971   6.890  1.00 41.07           C
ATOM    322  O   ARG H  48      80.596 -29.428   7.997  1.00 40.42           O
ATOM    323  CB  ARG H  48      80.994 -30.328   4.884  1.00 42.65           C
ATOM    324  CG  ARG H  48      81.469 -31.626   5.452  1.00 44.32           C
ATOM    325  CD  ARG H  48      81.043 -32.849   4.627  1.00 44.70           C
ATOM    326  NE  ARG H  48      81.529 -34.052   5.294  1.00 44.81           N
ATOM    327  CZ  ARG H  48      82.738 -34.580   5.140  1.00 45.11           C
ATOM    328  NH1 ARG H  48      83.614 -34.042   4.313  1.00 43.19           N
ATOM    329  NH2 ARG H  48      83.095 -35.620   5.877  1.00 47.83           N
ATOM    330  N   LEU H  49      79.221 -28.306   6.614  1.00 41.54           N
ATOM    331  CA  LEU H  49      78.167 -28.103   7.596  1.00 42.36           C
ATOM    332  C   LEU H  49      76.946 -28.837   7.064  1.00 41.91           C
ATOM    333  O   LEU H  49      76.756 -28.928   5.852  1.00 43.00           O
ATOM    334  CB  LEU H  49      77.839 -26.618   7.736  1.00 43.11           C
ATOM    335  CG  LEU H  49      78.845 -25.691   8.414  1.00 45.90           C
ATOM    336  CD1 LEU H  49      78.506 -24.251   8.063  1.00 46.63           C
ATOM    337  CD2 LEU H  49      78.809 -25.892   9.919  1.00 47.40           C
ATOM    338  N   GLY H  49A     76.120 -29.358   7.965  1.00 40.67           N
ATOM    339  CA  GLY H  49A     74.938 -30.081   7.530  1.00 39.48           C
ATOM    340  C   GLY H  49A     75.302 -31.190   6.558  1.00 39.18           C
ATOM    341  O   GLY H  49A     74.504 -31.580   5.700  1.00 38.72           O
ATOM    342  N   GLY H  50      76.527 -31.692   6.695  1.00 39.31           N
ATOM    343  CA  GLY H  50      77.002 -32.764   5.839  1.00 38.66           C
ATOM    344  C   GLY H  50      77.149 -32.388   4.377  1.00 38.02           C
ATOM    345  O   GLY H  50      77.152 -33.261   3.507  1.00 38.89           O
ATOM    346  N   ILE H  51      77.278 -31.095   4.097  1.00 36.09           N
ATOM    347  CA  ILE H  51      77.421 -30.635   2.725  1.00 33.47           C
ATOM    348  C   ILE H  51      78.689 -29.823   2.570  1.00 32.56           C
ATOM    349  O   ILE H  51      78.976 -28.949   3.378  1.00 31.71           O
ATOM    350  CB  ILE H  51      76.229 -29.776   2.323  1.00 32.60           C
ATOM    351  CG1 ILE H  51      74.948 -30.588   2.503  1.00 33.26           C
ATOM    352  CG2 ILE H  51      76.385 -29.303   0.880  1.00 31.47           C
ATOM    353  CD1 ILE H  51      73.684 -29.786   2.351  1.00 33.83           C
ATOM    354  N   ALA H  52      79.443 -30.112   1.520  1.00 31.45           N
ATOM    355  CA  ALA H  52      80.688 -29.412   1.265  1.00 31.24           C
ATOM    356  C   ALA H  52      80.477 -27.950   0.892  1.00 32.00           C
ATOM    357  O   ALA H  52      79.410 -27.556   0.406  1.00 32.53           O
ATOM    358  CB  ALA H  52      81.455 -30.109   0.148  1.00 27.90           C
'''

test_pdb_4 = '''\
CRYST1  577.812  448.715  468.790  90.00  90.00  90.00 P 1
SCALE1      0.001731  0.000000  0.000000        0.00000
SCALE2      0.000000  0.002229  0.000000        0.00000
SCALE3      0.000000  0.000000  0.002133        0.00000
ATOM      1  CA  LYS A 151      10.766   9.333  12.905  1.00 44.22           C
ATOM      2  CA  LYS A 152      10.117   9.159  11.610  1.00 49.42           C
ATOM      3  CA  LYS A 153       9.099   8.000  11.562  1.00 46.15           C
ATOM      4  CA  LYS A 154       8.000   8.202  11.065  1.00 52.97           C
ATOM      5  CA  LYS A 155      11.146   9.065  10.474  1.00 41.68           C
ATOM      6  CA  LYS A 156      10.547   9.007   9.084  1.00 55.55           C
ATOM      7  CA  LYS A 157      11.545   9.413   8.000  1.00 72.27           C
HETATM 3512  O   HOH A2001     -85.460  23.570 -79.414  1.00 34.05           O
HETATM 3513  O   HOH A2002     -81.492  18.186 -35.869  1.00 13.88           O
HETATM 3514  O   HOH A2003     -73.597  30.740 -72.170  1.00 33.13           O
HETATM 3515  O   HOH A2004     -74.933  23.929 -66.597  1.00 26.66           O
HETATM 3516  O   HOH A2005     -73.036  25.921 -66.605  1.00 23.41
END
'''

test_pdb_5 = '''\
CRYST1   42.558  113.820  137.132  90.00  90.00  90.00 P 21 21 21   24
SCALE1      0.023497  0.000000  0.000000        0.00000
SCALE2      0.000000  0.008786  0.000000        0.00000
SCALE3      0.000000  0.000000  0.007292        0.00000
ATOM   2905  N   ILE D  25     -21.903  62.095  51.650  1.00 19.95           N
ATOM   2906  CA  ILE D  25     -20.570  61.733  52.122  1.00 20.06           C
ATOM   2907  C   ILE D  25     -19.631  62.938  52.275  1.00 22.71           C
ATOM   2908  O   ILE D  25     -19.519  63.807  51.393  1.00 21.72           O
ATOM   2909  CB  ILE D  25     -19.887  60.655  51.241  1.00 15.34           C
ATOM   2910  CG1 ILE D  25     -20.768  59.411  51.256  1.00 12.56           C
ATOM   2911  CG2 ILE D  25     -18.438  60.365  51.831  1.00 14.73           C
ATOM   2912  CD1 ILE D  25     -20.661  58.435  50.096  1.00 15.09           C
ATOM   2913  N   GLY D  26     -18.968  62.975  53.429  1.00 20.40           N
ATOM   2914  CA  GLY D  26     -17.920  63.914  53.711  1.00 22.40           C
ATOM   2915  C   GLY D  26     -18.420  65.290  54.086  1.00 23.10           C
ATOM   2916  O   GLY D  26     -18.155  65.781  55.197  1.00 24.70           O
ATOM   2917  N   VAL D  27     -19.117  65.935  53.151  1.00 23.71           N
ATOM   2918  CA  VAL D  27     -19.595  67.307  53.388  1.00 23.17           C
ATOM   2919  C   VAL D  27     -20.597  67.344  54.568  1.00 24.78           C
ATOM   2920  O   VAL D  27     -20.745  68.387  55.199  1.00 26.51           O
ATOM   2921  CB  VAL D  27     -20.183  67.930  52.109  1.00 25.99           C
ATOM   2922  CG1 VAL D  27     -19.069  68.049  51.032  1.00 23.04           C
ATOM   2923  CG2 VAL D  27     -21.399  67.121  51.609  1.00 21.43           C
ATOM   2932  CB  ASP D  28     -23.555  65.658  55.450  1.00 25.04           C
ATOM   2933  CG  ASP D  28     -24.173  66.689  54.514  1.00 28.36           C
ATOM   2934  OD1 ASP D  28     -23.934  67.903  54.695  1.00 27.71           O
ATOM   2935  OD2 ASP D  28     -24.948  66.267  53.632  1.00 24.68           O
ATOM   2924  N  AASP D  28     -21.219  66.203  54.890  0.50 22.16           N
ATOM   2926  CA AASP D  28     -22.200  66.120  55.994  0.50 24.54           C
ATOM   2928  C  AASP D  28     -21.780  65.239  57.184  0.50 21.81           C
ATOM   2930  O  AASP D  28     -22.627  64.713  57.901  0.50 22.07           O
ATOM   2925  N  BASP D  28     -21.261  66.213  54.799  0.50 23.07           N
ATOM   2927  CA BASP D  28     -22.155  66.024  55.923  0.50 25.72           C
ATOM   2929  C  BASP D  28     -21.575  64.885  56.777  0.50 24.87           C
ATOM   2931  O  BASP D  28     -21.997  63.731  56.681  0.50 28.18           O
ATOM   2936  N   SER D  29     -20.479  65.193  57.444  1.00 26.01           N
ATOM   2937  CA  SER D  29     -19.892  64.371  58.495  1.00 24.17           C
ATOM   2938  C   SER D  29     -19.715  65.264  59.716  1.00 23.31           C
ATOM   2939  O   SER D  29     -19.972  66.479  59.622  1.00 23.57           O
ATOM   2940  CB  SER D  29     -18.538  63.871  58.037  1.00 21.39           C
ATOM   2941  OG  SER D  29     -17.658  64.927  57.738  1.00 20.49           O
ATOM   2942  N   ALA D  30     -19.259  64.708  60.836  1.00 19.54           N
ATOM   2943  CA  ALA D  30     -18.937  65.524  61.993  1.00 21.15           C
ATOM   2944  C   ALA D  30     -17.875  66.586  61.555  1.00 20.92           C
ATOM   2945  O   ALA D  30     -18.026  67.767  61.866  1.00 18.58           O
ATOM   2946  CB  ALA D  30     -18.461  64.657  63.163  1.00 17.17           C
'''

test_pdb_6 = '''\
CRYST1  203.106   83.279  178.234  90.00 106.67  90.00 C 1 2 1      12
ATOM      1  N   ARG H  48      82.680 -29.174   6.224  1.00 40.63           N
ATOM      2  CA  ARG H  48      81.318 -29.102   5.747  1.00 40.99           C
ATOM      3  C   ARG H  48      80.335 -28.971   6.890  1.00 41.07           C
ATOM      4  O   ARG H  48      80.596 -29.428   7.997  1.00 40.42           O
ATOM      5  N   LEU H  49      79.221 -28.306   6.614  1.00 41.54           N
ATOM      6  CA  LEU H  49      78.167 -28.103   7.596  1.00 42.36           C
ATOM      7  C   LEU H  49      76.946 -28.837   7.064  1.00 41.91           C
ATOM      8  O   LEU H  49      76.756 -28.928   5.852  1.00 43.00           O
ATOM      9  CB  LEU H  49      77.839 -26.618   7.736  1.00 43.11           C
ATOM     10  CG  LEU H  49      78.845 -25.691   8.414  1.00 45.90           C
ATOM     11  CD1 LEU H  49      78.506 -24.251   8.063  1.00 46.63           C
ATOM     12  CD2 LEU H  49      78.809 -25.892   9.919  1.00 47.40           C
ATOM     13  N   GLY H  49A     76.120 -29.358   7.965  1.00 40.67           N
ATOM     14  CA  GLY H  49A     74.938 -30.081   7.530  1.00 39.48           C
ATOM     15  C   GLY H  49A     75.302 -31.190   6.558  1.00 39.18           C
ATOM     16  O   GLY H  49A     74.504 -31.580   5.700  1.00 38.72           O
ATOM     17  N   GLY H  50      76.527 -31.692   6.695  1.00 39.31           N
ATOM     18  CA  GLY H  50      77.002 -32.764   5.839  1.00 38.66           C
ATOM     19  C   GLY H  50      77.149 -32.388   4.377  1.00 38.02           C
ATOM     20  O   GLY H  50      77.152 -33.261   3.507  1.00 38.89           O
'''

test_pdb_7 = """\
ATOM      1  N   ILE A  14     -14.366 134.691 190.492  1.00 65.50           N
ATOM      2  CA  ILE A  14     -13.578 133.563 190.974  1.00 63.49           C
ATOM      3  C   ILE A  14     -12.354 134.128 191.678  1.00 62.68           C
ATOM      4  O   ILE A  14     -11.609 134.922 191.093  1.00 63.22           O
ATOM      5  CB  ILE A  14     -13.151 132.610 189.839  1.00 63.05           C
ATOM      6  CG1 ILE A  14     -14.374 131.920 189.230  1.00 63.83           C
ATOM      7  CG2 ILE A  14     -12.167 131.562 190.366  1.00 61.34           C
ATOM      8  CD1 ILE A  14     -14.061 131.063 188.013  1.00 63.78           C
ATOM      9  N   ALA A  15     -12.147 133.721 192.927  1.00 60.93           N
ATOM     10  CA  ALA A  15     -11.017 134.168 193.727  1.00 60.26           C
ATOM     11  C   ALA A  15     -10.262 132.961 194.262  1.00 58.80           C
ATOM     12  O   ALA A  15     -10.871 131.958 194.644  1.00 58.33           O
ATOM     13  CB  ALA A  15     -11.474 135.043 194.893  1.00 60.84           C
ATOM     14  N   GLY A  16      -8.937 133.066 194.291  1.00 59.81           N
ATOM     15  CA  GLY A  16      -8.107 131.990 194.797  1.00 58.83           C
ATOM     16  C   GLY A  16      -6.893 132.537 195.514  1.00 58.70           C
ATOM     17  O   GLY A  16      -6.463 133.670 195.277  1.00 59.33           O
ATOM     18  N   ASP A  17      -6.340 131.714 196.401  1.00 58.04           N
ATOM     19  CA  ASP A  17      -5.148 132.100 197.141  1.00 58.08           C
ATOM     20  C   ASP A  17      -3.917 131.803 196.282  1.00 58.37           C
ATOM     21  O   ASP A  17      -4.024 131.476 195.097  1.00 58.58           O
ATOM     22  CB  ASP A  17      -5.110 131.397 198.496  1.00 57.69           C
ATOM     23  CG  ASP A  17      -4.577 129.986 198.405  1.00 57.48           C
ATOM     24  OD1 ASP A  17      -4.995 129.243 197.491  1.00 57.38           O
ATOM     25  OD2 ASP A  17      -3.733 129.619 199.249  1.00 57.63           O1-
"""

test_pdb_8 = """\
ATOM      1  N   ALA A 125     -29.026  54.443 -42.502  1.00 32.57           N
ATOM      2  CA  ALA A 125     -27.931  53.547 -42.137  1.00 31.14           C
ATOM      3  C   ALA A 125     -28.182  52.145 -42.679  1.00 31.38           C
ATOM      4  O   ALA A 125     -27.272  51.499 -43.215  1.00 39.47           O
ATOM      5  CB  ALA A 125     -27.774  53.497 -40.615  1.00 27.41           C
ATOM      6  N   SER A 126     -29.423  51.670 -42.564  1.00 30.06           N
ATOM      7  CA  SER A 126     -29.767  50.337 -43.032  1.00 32.18           C
ATOM      8  C   SER A 126     -29.602  50.216 -44.538  1.00 38.36           C
ATOM      9  O   SER A 126     -29.386  49.109 -45.039  1.00 45.79           O
ATOM     10  CB  SER A 126     -31.190  49.988 -42.605  1.00 29.70           C
ATOM     11  OG  SER A 126     -31.370  48.579 -42.545  1.00 40.98           O
ATOM     12  N   ALA A 127     -29.671  51.312 -45.252  1.00 45.06           N
ATOM     13  CA  ALA A 127     -29.521  51.369 -46.685  1.00 48.16           C
ATOM     14  C   ALA A 127     -28.142  51.136 -47.248  1.00 53.02           C
ATOM     15  O   ALA A 127     -28.024  50.717 -48.363  1.00 58.10           O
ATOM     16  CB  ALA A 127     -30.098  52.647 -47.231  1.00 47.36           C
ATOM     17  N   LYS A 128     -27.105  51.433 -46.492  1.00 54.75           N
ATOM     18  CA  LYS A 128     -25.742  51.247 -46.964  1.00 57.26           C
ATOM     19  C   LYS A 128     -25.482  49.795 -47.340  1.00 52.16           C
ATOM     20  O   LYS A 128     -25.996  48.893 -46.720  1.00 45.54           O
ATOM     21  CB  LYS A 128     -24.759  51.535 -45.843  1.00 30.00           C
ATOM     22  CG  LYS A 128     -24.496  52.960 -45.461  1.00 30.00           C
ATOM     23  CD  LYS A 128     -22.989  53.226 -45.366  1.00 30.00           C
ATOM     24  CE  LYS A 128     -22.248  52.508 -44.239  1.00 30.00           C
ATOM     25  NZ  LYS A 128     -20.903  53.077 -43.987  1.00 30.00           N
"""

test_pdb_9 = """\
ATOM      1  N   ARG A 124     -31.625  55.204 -43.320  1.00 36.06           N
ATOM      2  CA  ARG A 124     -30.475  55.601 -44.092  1.00 34.32           C
ATOM      3  C   ARG A 124     -29.278  54.709 -43.767  1.00 36.61           C
ATOM      4  O   ARG A 124     -28.616  54.266 -44.645  1.00 40.69           O
ATOM      5  CB  ARG A 124     -30.122  57.063 -43.899  1.00 30.00           C
ATOM      6  CG  ARG A 124     -28.740  57.425 -44.418  1.00 30.00           C
ATOM      7  CD  ARG A 124     -28.283  58.842 -44.132  1.00 30.00           C
ATOM      8  NE  ARG A 124     -28.977  59.388 -42.992  1.00 30.00           N
ATOM      9  CZ  ARG A 124     -28.813  60.609 -42.544  1.00 30.00           C
ATOM     10  NH1 ARG A 124     -27.970  61.414 -43.148  1.00 30.00           N
ATOM     11  NH2 ARG A 124     -29.492  61.020 -41.504  1.00 30.00           N
ATOM     12  N   ALA A 125     -29.026  54.443 -42.502  1.00 32.57           N
ATOM     13  CA  ALA A 125     -27.931  53.547 -42.137  1.00 31.14           C
ATOM     14  C   ALA A 125     -28.182  52.145 -42.679  1.00 31.38           C
ATOM     15  O   ALA A 125     -27.272  51.499 -43.215  1.00 39.47           O
ATOM     16  CB  ALA A 125     -27.774  53.497 -40.615  1.00 27.41           C
ATOM     17  N   SER A 126     -29.423  51.670 -42.564  1.00 30.06           N
ATOM     18  CA  SER A 126     -29.767  50.337 -43.032  1.00 32.18           C
ATOM     19  C   SER A 126     -29.602  50.216 -44.538  1.00 38.36           C
ATOM     20  O   SER A 126     -29.386  49.109 -45.039  1.00 45.79           O
ATOM     21  CB  SER A 126     -31.190  49.988 -42.605  1.00 29.70           C
ATOM     22  OG  SER A 126     -31.370  48.579 -42.545  1.00 40.98           O
ATOM     23  N   ALA A 127     -29.671  51.312 -45.252  1.00 45.06           N
ATOM     24  CA  ALA A 127     -29.521  51.369 -46.685  1.00 48.16           C
ATOM     25  C   ALA A 127     -28.142  51.136 -47.248  1.00 53.02           C
ATOM     26  O   ALA A 127     -28.024  50.717 -48.363  1.00 58.10           O
ATOM     27  CB  ALA A 127     -30.098  52.647 -47.231  1.00 47.36           C
"""

test_pdb_10 = """\
ATOM    885  N   ASP A 117     -38.120  61.598 -35.985  1.00 24.56           N
ATOM    886  CA  ASP A 117     -37.446  61.205 -37.186  1.00 24.60           C
ATOM    887  C   ASP A 117     -36.374  60.119 -37.011  1.00 26.59           C
ATOM    888  O   ASP A 117     -35.562  59.977 -37.870  1.00 19.99           O
ATOM    889  CB  ASP A 117     -36.851  62.414 -37.878  1.00 33.00           C
ATOM    890  CG  ASP A 117     -37.507  62.725 -39.182  1.00 46.28           C
ATOM    891  OD1 ASP A 117     -38.553  62.190 -39.479  1.00 53.20           O
ATOM    892  OD2 ASP A 117     -36.964  63.519 -39.921  1.00 54.16           O
ATOM      1  N   ARG A 124     -31.625  55.204 -43.320  1.00 36.06           N
ATOM      2  CA  ARG A 124     -30.475  55.601 -44.092  1.00 34.32           C
ATOM      3  C   ARG A 124     -29.278  54.709 -43.767  1.00 36.61           C
ATOM      4  O   ARG A 124     -28.616  54.266 -44.645  1.00 40.69           O
ATOM      5  CB  ARG A 124     -30.122  57.063 -43.899  1.00 30.00           C
ATOM      6  CG  ARG A 124     -28.740  57.425 -44.418  1.00 30.00           C
ATOM      7  CD  ARG A 124     -28.283  58.842 -44.132  1.00 30.00           C
ATOM      8  NE  ARG A 124     -28.977  59.388 -42.992  1.00 30.00           N
ATOM      9  CZ  ARG A 124     -28.813  60.609 -42.544  1.00 30.00           C
ATOM     10  NH1 ARG A 124     -27.970  61.414 -43.148  1.00 30.00           N
ATOM     11  NH2 ARG A 124     -29.492  61.020 -41.504  1.00 30.00           N
ATOM     12  N   ALA A 125     -29.026  54.443 -42.502  1.00 32.57           N
ATOM     13  CA  ALA A 125     -27.931  53.547 -42.137  1.00 31.14           C
ATOM     14  C   ALA A 125     -28.182  52.145 -42.679  1.00 31.38           C
ATOM     15  O   ALA A 125     -27.272  51.499 -43.215  1.00 39.47           O
ATOM     16  CB  ALA A 125     -27.774  53.497 -40.615  1.00 27.41           C
ATOM     17  N   SER A 126     -29.423  51.670 -42.564  1.00 30.06           N
ATOM     18  CA  SER A 126     -29.767  50.337 -43.032  1.00 32.18           C
ATOM     19  C   SER A 126     -29.602  50.216 -44.538  1.00 38.36           C
ATOM     20  O   SER A 126     -29.386  49.109 -45.039  1.00 45.79           O
ATOM     21  CB  SER A 126     -31.190  49.988 -42.605  1.00 29.70           C
ATOM     22  OG  SER A 126     -31.370  48.579 -42.545  1.00 40.98           O
ATOM     23  N   ALA A 127     -29.671  51.312 -45.252  1.00 45.06           N
ATOM     24  CA  ALA A 127     -29.521  51.369 -46.685  1.00 48.16           C
ATOM     25  C   ALA A 127     -28.142  51.136 -47.248  1.00 53.02           C
ATOM     26  O   ALA A 127     -28.024  50.717 -48.363  1.00 58.10           O
ATOM     27  CB  ALA A 127     -30.098  52.647 -47.231  1.00 47.36           C
"""

test_pdb_11 = """\
ATOM      1  CA  MET A 480     -10.152  -1.677  37.457  1.00 28.90      AA-  C
ATOM      2  CG  MET A 480      -9.517   0.607  38.489  1.00 49.21      AA-  C
ATOM      3  CA  VAL A 481     -11.857  -5.041  37.936  1.00 23.15      AA-  C
ATOM      4  CA  VAL A 482     -13.580  -5.715  41.251  1.00 34.23      AA-  C
ATOM      5  CA  ARG A 483     -14.378  -9.189  42.566  1.00 30.50      AA-  C
ATOM      6  CG  ARG A 483     -12.194 -10.275  42.051  1.00 28.52      AA-  C
ATOM      7  CA  ARG A 484     -17.677  -9.381  44.419  1.00 41.14      AA-  C
ATOM      8  CA  MET B 480     -16.603   9.979   4.769  1.00 32.68      BA-  C
ATOM      9  CG  MET B 480     -16.554   7.723   3.530  1.00 54.69      BA-  C
ATOM     10  CA  VAL B 481     -17.675  13.614   4.970  1.00 34.88      BA-  C
ATOM     11  CA  VAL B 482     -20.607  14.715   2.834  1.00 35.95      BA-  C
ATOM     12  CA  ARG B 483     -21.481  18.297   1.835  1.00 23.95      BA-  C
ATOM     13  CA  ARG B 484     -25.162  19.109   1.469  1.00 47.03      BA-  C
ATOM     14  CG  ARG B 484     -25.746  17.004   2.827  1.00 47.49      BA-  C
"""

test_pdb_12 = """\
ATOM      1  N   ASP A 279     148.294 135.256  14.514  1.00361.90           N
ATOM      2  CA  ASP A 279     148.944 135.688  13.244  1.00363.91           C
ATOM      3  C   ASP A 279     150.116 136.660  13.498  1.00311.37           C
ATOM      4  O   ASP A 279     150.890 136.977  12.584  1.00265.25           O
ATOM      5  CB  ASP A 279     147.909 136.283  12.247  1.00340.39           C
ATOM      6  CG  ASP A 279     148.544 136.866  10.940  1.00317.09           C
ATOM      7  OD1 ASP A 279     149.462 136.253  10.337  1.00271.02           O
ATOM      8  OD2 ASP A 279     148.085 137.945  10.494  1.00309.85           O
ATOM      9  N   ARG A 280     150.240 137.122  14.743  1.00282.66           N
ATOM     10  CA  ARG A 280     151.405 137.889  15.191  1.00243.53           C
ATOM     11  C   ARG A 280     152.420 136.982  15.887  1.00239.61           C
ATOM     12  O   ARG A 280     153.345 137.464  16.539  1.00231.28           O
ATOM     13  CB  ARG A 280     150.976 139.011  16.137  1.00252.04           C
ATOM     14  N   LEU A 281     152.236 135.669  15.740  1.00271.48           N
ATOM     15  CA  LEU A 281     153.120 134.676  16.348  1.00298.95           C
ATOM     16  C   LEU A 281     153.688 133.669  15.334  1.00293.08           C
ATOM     17  O   LEU A 281     154.558 132.864  15.675  1.00325.01           O
ATOM     18  CB  LEU A 281     152.411 133.959  17.495  1.00294.06           C
ATOM     19  CG  LEU A 281     153.236 133.921  18.776  1.00316.57           C
ATOM     20  CD1 LEU A 281     153.263 135.296  19.436  1.00243.06           C
ATOM     21  CD2 LEU A 281     152.663 132.873  19.709  1.00298.46           C
ATOM     22  N   GLU A 282     153.167 133.706  14.106  1.00247.81           N
ATOM     23  CA  GLU A 282     153.827 133.132  12.921  1.00294.79           C
ATOM     24  C   GLU A 282     154.887 134.153  12.467  1.00320.33           C
ATOM     25  O   GLU A 282     155.953 133.788  11.917  1.00233.13           O
ATOM     26  CB  GLU A 282     152.782 132.842  11.808  1.00321.01           C
ATOM     27  CG  GLU A 282     153.279 132.850  10.343  1.00432.99           C
ATOM     28  CD  GLU A 282     152.275 133.402   9.306  1.00369.39           C
ATOM     29  OE1 GLU A 282     151.087 133.604   9.642  1.00323.25           O
ATOM     30  OE2 GLU A 282     152.674 133.646   8.135  1.00284.42           O
ATOM     31  N   ARG A 283     154.572 135.429  12.728  1.00357.16           N
ATOM     32  CA  ARG A 283     155.454 136.572  12.483  1.00384.53           C
ATOM     33  C   ARG A 283     156.436 136.825  13.640  1.00275.79           C
ATOM     34  O   ARG A 283     157.288 137.710  13.545  1.00265.14           O
ATOM     35  CB  ARG A 283     154.626 137.840  12.208  1.00316.05           C
ATOM     36  N   ARG A 284     156.321 136.051  14.722  1.00266.18           N
ATOM     37  CA  ARG A 284     157.243 136.156  15.872  1.00255.82           C
ATOM     38  C   ARG A 284     158.154 134.935  16.094  1.00314.62           C
ATOM     39  O   ARG A 284     159.154 135.024  16.825  1.00253.48           O
ATOM     40  CB  ARG A 284     156.469 136.452  17.153  1.00230.14           C
ATOM     41  N   SER A 285     157.800 133.808  15.471  1.00393.24           N
ATOM     42  CA  SER A 285     158.572 132.562  15.569  1.00358.40           C
ATOM     43  C   SER A 285     159.419 132.307  14.317  1.00418.23           C
ATOM     44  O   SER A 285     158.898 132.157  13.204  1.00380.17           O
ATOM     45  CB  SER A 285     157.651 131.363  15.841  1.00264.20           C
"""

test_pdb_13 = """\
ATOM      1  N   GLU A 260      74.225  41.166  25.753  1.00 58.41           N
ATOM      2  CA  GLU A 260      75.147  41.741  24.757  1.00 58.93           C
ATOM      3  C   GLU A 260      75.096  41.010  23.404  1.00 59.55           C
ATOM      4  O   GLU A 260      75.896  41.322  22.519  1.00 59.69           O
ATOM      5  CB  GLU A 260      76.587  41.677  25.279  1.00 59.13           C
ATOM      6  CG  GLU A 260      76.899  42.547  26.490  1.00 59.98           C
ATOM      7  CD  GLU A 260      78.326  42.347  27.001  1.00 61.43           C
ATOM      8  OE1 GLU A 260      78.859  41.213  26.938  1.00 61.60           O
ATOM      9  OE2 GLU A 260      78.934  43.331  27.468  1.00 64.79           O
ATOM     10  N   CYS A 261      74.177  40.043  23.246  1.00 60.12           N
ATOM     11  CA  CYS A 261      74.105  39.192  22.050  1.00 60.18           C
ATOM     12  C   CYS A 261      72.733  39.239  21.359  1.00 61.64           C
ATOM     13  O   CYS A 261      72.355  38.276  20.698  1.00 61.42           O
ATOM     14  CB  CYS A 261      74.439  37.718  22.407  1.00 59.87           C
ATOM     15  SG  CYS A 261      76.011  37.403  23.356  1.00 57.71           S
ATOM     16  N   GLY A 262      72.018  40.364  21.494  1.00 63.62           N
ATOM     17  CA  GLY A 262      70.666  40.587  20.921  1.00 65.14           C
ATOM     18  C   GLY A 262      69.546  40.080  21.841  1.00 66.01           C
ATOM     19  O   GLY A 262      68.528  39.491  21.408  1.00 66.75           O
ATOM     20  OXT GLY A 262      69.622  40.268  23.068  1.00 66.78           O
TER
ATOM     21  N   CYS B 261      74.710  74.423  33.026  1.00 51.47           N
ATOM     22  CA  CYS B 261      74.050  74.094  31.771  1.00 51.64           C
ATOM     23  C   CYS B 261      73.663  75.409  31.099  1.00 53.38           C
ATOM     24  O   CYS B 261      72.609  75.498  30.502  1.00 53.30           O
ATOM     25  CB  CYS B 261      72.820  73.215  32.017  1.00 50.99           C
ATOM     26  SG  CYS B 261      73.155  71.752  33.046  1.00 46.82           S
ATOM     27  N   GLY B 262      74.546  76.413  31.226  1.00 55.76           N
ATOM     28  CA  GLY B 262      74.424  77.787  30.642  1.00 57.82           C
ATOM     29  C   GLY B 262      73.715  78.804  31.551  1.00 58.64           C
ATOM     30  O   GLY B 262      73.579  78.604  32.783  1.00 58.58           O
ATOM     31  OXT GLY B 262      73.311  79.889  31.060  1.00 59.70           O
TER
HETATM   32  N   GLU A 301      75.002  26.783  43.217  1.00 28.86           N
HETATM   33  CA  GLU A 301      76.076  25.770  43.503  1.00 29.91           C
HETATM   34  C   GLU A 301      76.194  25.508  45.002  1.00 30.65           C
HETATM   35  O   GLU A 301      77.025  24.681  45.464  1.00 30.01           O
HETATM   36  CB  GLU A 301      77.422  26.205  42.939  1.00 30.22           C
HETATM   37  CG  GLU A 301      77.563  26.074  41.392  1.00 29.97           C
HETATM   38  CD  GLU A 301      77.627  24.641  40.897  1.00 30.03           C
HETATM   39  OE1 GLU A 301      77.884  23.761  41.746  1.00 26.02           O
HETATM   40  OE2 GLU A 301      77.475  24.404  39.651  1.00 30.03           O
HETATM   41  OXT GLU A 301      75.416  26.080  45.793  1.00 30.21           O
"""

def test_get_clean_selection_string():
  """ Check get_clean_selection_string  """
  # print sys._getframe().f_code.co_name
  ch_sel = 'chain A'
  res_selection1 = ['resid 1:10']

  res_id = '27c'
  s = '(resid ' + res_id + ' and (name '
  atom_name = [x for x in [' CA',' N']]
  atom_str = ' or name '.join(atom_name)
  s = s + atom_str + ')'
  res_selection2 = ['resid 1:10',s]

  result1 = get_clean_selection_string(ch_sel,[])
  result2 = get_clean_selection_string(ch_sel,res_selection1)
  result3 = get_clean_selection_string(ch_sel,res_selection2)

  assert result1 == 'chain A', result1
  assert result2 == 'chain A and resid 1:10', result2
  expt = 'chain A and (resid 1:10 or (resid 27c and (name CA or name N))'
  assert result3 == expt, result3

def test_selection_string_from_selection():
  """ Test selection_string_from_selection """
  pdb_h = iotbx.pdb.input(source_info=None, lines=test_pdb_1).construct_hierarchy()
  isel1 = flex.size_t([12, 13, 14, 15, 16, 17, 18])
  isel2 = flex.size_t([12, 13, 14, 16, 17, 18])
  isel3 = flex.size_t([12, 13, 14, 15, 16, 17])
  sel_str1 = selection_string_from_selection(pdb_h,isel1)
  sel_str2 = selection_string_from_selection(pdb_h,isel2)
  sel_str3 = selection_string_from_selection(pdb_h,isel3)
  assert sel_str1 == "chain 'D'", sel_str1
  assert sel_str2 == "(chain 'D' and (resid 1 through 3 or resid 5 through 7))", sel_str2
  assert sel_str3 == "(chain 'D' and resid 1 through 6)", sel_str3
  #
  asc = pdb_h.atom_selection_cache()
  sel1 = list(asc.iselection(sel_str1))
  sel2 = list(asc.iselection(sel_str2))
  sel3 = list(asc.iselection(sel_str3))
  #
  assert sel1 == list(isel1), sel1
  assert sel2 == list(isel2), sel2
  assert sel3 == list(isel3), sel3

def test_selection_string_from_selection2():
  """ Test selection_string_from_selection """
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_2).construct_hierarchy()
  l1 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
        26,27,28,29,30,31]
  l2 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17      ,20,21,22,23,24,25,
        26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42]
  isel1 = flex.size_t(l1)
  isel2 = flex.size_t(l2)
  #
  sel_str1 = selection_string_from_selection(pdb_h,isel1)
  sel_str2 = selection_string_from_selection(pdb_h,isel2)
  assert sel_str1 == "chain 'A' or (chain 'B' and resid 153 through 154)", sel_str1
  s = "(chain 'A' and (resid 151 or (resid 152 and (name N or name CA or "
  s += "name C or name O or name CB or name CG or name CD or name NE or "
  s += "name CZ )))) or chain 'B'"
  assert sel_str2 == s, sel_str2
  #
  asc = pdb_h.atom_selection_cache()
  sel1 = list(asc.iselection(sel_str1))
  sel2 = list(asc.iselection(sel_str2))
  #
  assert sel1 == list(isel1), sel1
  assert sel2 == list(isel2), sel2

def test_avoid_chain_selection():
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_2).construct_hierarchy()
  isel1 = flex.size_t([0,1,2,3,4,5,6,7,8])
  sel_str1 = selection_string_from_selection(pdb_h,isel1)
  s = "(chain 'A' and resid 151)"
  assert sel_str1 == s, sel_str1

def test_avoid_chain_selection2():
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_3).construct_hierarchy()
  isel1 = flex.size_t(range(6,46))
  sel_str1 = selection_string_from_selection(pdb_h,isel1)
  # s = '(chain H and (resid 48 or resid 49 or resid 49A or resid 50:52))'
  # better way:
  s = "(chain 'H' and resid 48 through 52)"
  assert sel_str1 == s, sel_str1
  #
  l1 = list(range(6,25)) + list(range(29,46))
  isel1 = flex.size_t(l1)
  # s = '(chain H and (resid 48 or resid 49 or resid 50:52))'
  # better way:
  s = "(chain 'H' and (resid 48 through 49 or resid 50 through 52))"
  sel_str1 = selection_string_from_selection(pdb_h,isel1)
  assert sel_str1 == s, sel_str1

def test_avoid_hoh():
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_4).construct_hierarchy()
  isel1 = flex.size_t(range(7))
  sel_str1 = selection_string_from_selection(pdb_h,isel1)
  s = "(chain 'A' and resid 151 through 157)"
  assert sel_str1 == s, sel_str1
  #
  asc = pdb_h.atom_selection_cache()
  sel = asc.iselection(s)
  assert sel.size() == 7, sel.size()

def test_include_hoh():
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_4).construct_hierarchy()
  isel1 = flex.size_t(range(7))
  sel_str1 = selection_string_from_selection(
    pdb_h,isel1)
  s = "(chain 'A' and resid 151 through 157)"
  assert sel_str1 == s, sel_str1
  #
  asc = pdb_h.atom_selection_cache()
  sel = asc.iselection(s)
  assert sel.size() == 7, sel.size()
  #
  isel1 = flex.size_t(range(12))
  sel_str1 = selection_string_from_selection(
    pdb_h,isel1)
  assert sel_str1 == "chain 'A'", sel_str

def test_selection_with_alternative_conformers():
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_5).construct_hierarchy(sort_atoms=True)
  asc = pdb_h.atom_selection_cache()
  chains_info = get_chains_info(pdb_h)
  ch_D = chains_info['D']
  # test conditions verification
  assert ch_D.no_altloc == [True, True, True, False, True, True]
  select_all = sorted([x for xi in ch_D.atom_selection for x in xi])
  test_list = list(asc.iselection("not altloc B"))
  assert select_all == test_list, "%s" % select_all

def test_insertions():
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_6).construct_hierarchy()
  isel = flex.size_t(range(15))
  tsel = selection_string_from_selection(pdb_h, isel)
  assert tsel == "(chain 'H' and (resid 48 through 49 or (resid 49A and (name N or name CA or name C ))))", tsel

  isel = flex.size_t(range(16))
  tsel = selection_string_from_selection(pdb_h, isel)
  assert tsel == "(chain 'H' and resid 48 through 49A)", tsel

def test_2():
  """
  behavior with GLY: don't stop selection string:
  if a user wants only N, CA, C, O atoms, there is no reason to break the
  selection range just because there is GLY and it doesn't need names.
  And we even skip the residue range here, tested extensively in test_5
  """
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_7).construct_hierarchy()
  isel = flex.size_t([0,1,2,3,8,9,10,11,13,14,15,16,17,18,19,20])
  tsel = selection_string_from_selection(pdb_h, isel)
  assert tsel == "(chain 'A' and (name N or name CA or name C or name O ))" , tsel

def test_3():
  """
  single atom selections
  """
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_6).construct_hierarchy()
  for i, answ in zip(range(20), [
      "(chain 'H' and (resid 48 and (name N )))",
      "(chain 'H' and (resid 48 and (name CA )))",
      "(chain 'H' and (resid 48 and (name C )))",
      "(chain 'H' and (resid 48 and (name O )))",
      "(chain 'H' and (resid 49 and (name N )))",
      "(chain 'H' and (resid 49 and (name CA )))",
      "(chain 'H' and (resid 49 and (name C )))",
      "(chain 'H' and (resid 49 and (name O )))",
      "(chain 'H' and (resid 49 and (name CB )))",
      "(chain 'H' and (resid 49 and (name CG )))",
      "(chain 'H' and (resid 49 and (name CD1)))",
      "(chain 'H' and (resid 49 and (name CD2)))",
      "(chain 'H' and (resid 49A and (name N )))",
      "(chain 'H' and (resid 49A and (name CA )))",
      "(chain 'H' and (resid 49A and (name C )))",
      "(chain 'H' and (resid 49A and (name O )))",
      "(chain 'H' and (resid 50 and (name N )))",
      "(chain 'H' and (resid 50 and (name CA )))",
      "(chain 'H' and (resid 50 and (name C )))",
      "(chain 'H' and (resid 50 and (name O )))"]):
    isel = flex.size_t([i])
    tsel = selection_string_from_selection(pdb_h, isel)
    assert tsel == answ, "%s != %s" % (tsel, answ)

def test_4():
  """
  double atoms selections
  """
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_6).construct_hierarchy()
  for i, answ in zip(range(0,20,2), [
      "(chain 'H' and (resid 48 and (name N or name CA )))",
      "(chain 'H' and (resid 48 and (name C or name O )))",
      "(chain 'H' and (resid 49 and (name N or name CA )))",
      "(chain 'H' and (resid 49 and (name C or name O )))",
      "(chain 'H' and (resid 49 and (name CB or name CG )))",
      "(chain 'H' and (resid 49 and (name CD1 or name CD2)))",
      "(chain 'H' and (resid 49A and (name N or name CA )))",
      "(chain 'H' and (resid 49A and (name C or name O )))",
      "(chain 'H' and (resid 50 and (name N or name CA )))",
      "(chain 'H' and (resid 50 and (name C or name O )))"]):
    isel = flex.size_t([i,i+1])
    tsel = selection_string_from_selection(pdb_h, isel)
    assert tsel == answ, "%s != %s" % (tsel, answ)
  # and now odd:
  for i, answ in zip(range(1,19,2), [
      "(chain 'H' and (resid 48 and (name CA or name C )))",
      "(chain 'H' and ((resid 48 and (name O )) or (resid 49 and (name N ))))",
      "(chain 'H' and (resid 49 and (name CA or name C )))",
      "(chain 'H' and (resid 49 and (name O or name CB )))",
      "(chain 'H' and (resid 49 and (name CG or name CD1)))",
      "(chain 'H' and ((resid 49 and (name CD2)) or (resid 49A and (name N ))))",
      "(chain 'H' and (resid 49A and (name CA or name C )))",
      "(chain 'H' and ((resid 49A and (name O )) or (resid 50 and (name N ))))",
      "(chain 'H' and (resid 50 and (name CA or name C )))"]):
    isel = flex.size_t([i,i+1])
    tsel = selection_string_from_selection(pdb_h, isel)
    assert tsel == answ, "%s != %s" % (tsel, answ)

def test_5():
  """
  Don't output the residue range if it covers whole chain, even if
  there is atom name selection
  """
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_7).construct_hierarchy()
  isel = flex.size_t(range(25))
  tsel = selection_string_from_selection(pdb_h, isel)
  assert tsel == "chain 'A'", tsel
  isel = flex.size_t([0,8,13,17])
  tsel = selection_string_from_selection(pdb_h, isel)
  assert tsel == "(chain 'A' and (name N ))", tsel
  isel = flex.size_t([0,1,8,9,13,14,17,18])
  tsel = selection_string_from_selection(pdb_h, isel)
  assert tsel == "(chain 'A' and (name N or name CA ))", tsel
  isel = flex.size_t([0,1,2,8,9,13,14,17,18])
  tsel = selection_string_from_selection(pdb_h, isel)
  assert tsel == "(chain 'A' and ((resid 14 and (name N or name CA or name C )) or (resid 15 through 17 and (name N or name CA ))))", tsel

def test_6():
  """
  previous range is all atoms selected, next residue is not, but selected
  atoms are the same as for the last residue. In this case the range with all
  atoms should be dumped.
  """
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_8).construct_hierarchy()
  isel = flex.size_t(range(21))
  tsel = selection_string_from_selection(pdb_h, isel)
  # print "tsel", tsel
  assert tsel == "(chain 'A' and (resid 125 through 127 or (resid 128 and (name N or name CA or name C or name O or name CB ))))" , tsel

def test_7():
  """
  """
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_9).construct_hierarchy()
  isel = flex.size_t([0,1,2,3,4]+list(range(11,27)))
  tsel = selection_string_from_selection(pdb_h, isel)
  assert tsel == "(chain 'A' and ((resid 124 through 125 and (name N or name CA or name C or name O or name CB )) or resid 126 through 127))", tsel
  # print "tsel", tsel

def test_8():
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_10).construct_hierarchy()
  isel = flex.size_t(list(range(8))+[8,9,10,11,12]+list(range(19,35)))
  tsel = selection_string_from_selection(pdb_h, isel)
  assert tsel == "(chain 'A' and (resid 117 or (resid 124 through 125 and (name N or name CA or name C or name O or name CB )) or resid 126 through 127))", tsel


def test_11():
  """
  outputting name selection at the end of hierarchy?..
  """
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_11).construct_hierarchy()
  isel = flex.size_t(list(range(5))+[6])
  tsel = selection_string_from_selection(pdb_h, isel)
  assert tsel == "(chain 'A' and (resid 480 through 482 or (resid 483 through 484 and (name CA ))))", tsel

def test_12():
  """
  Not the first range in hierarchy, some atoms are absent for the first residue
  in the range, but atoms of the several next residues are coniside with
  present atoms of the first residue. And hierarchy ends. Make sure list
  of atoms is outputted for the last range.
  """
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_12).construct_hierarchy()
  isel = flex.size_t(list(range(26))+list(range(30,45)))
  tsel = selection_string_from_selection(pdb_h, isel)
  assert tsel == "(chain 'A' and (resid 279 through 281 or (resid 282 through 285 and (name N or name CA or name C or name O or name CB ))))", tsel

def test_13():
  """
  Test correct handling of portion of a chain at the end of the pdb file.
  """
  pdb_h = iotbx.pdb.input(
      source_info=None, lines=test_pdb_13).construct_hierarchy()
  isel = flex.size_t(list(range(20))+list(range(31,41)))
  tsel = selection_string_from_selection(pdb_h, isel)
  # print tsel
  assert tsel == "chain 'A'"
  isel = flex.size_t(list(range(19))+list(range(31,41)))
  tsel = selection_string_from_selection(pdb_h, isel)
  # print tsel
  assert tsel == "(chain 'A' and (resid 260 through 261 or (resid 262 and (name N or name CA or name C or name O )) or resid 301))"

def test_14():
  """
  Test correct handling of altloc upper and lower cases.
  """
  lines = """
CRYST1   45.640   40.754   30.275  90.00  90.00  90.00 P 21 21 21
ATOM      2  CA AVAL A   1      -4.890   1.653  12.849  0.50 13.49           C
ATOM      9  CA aVAL A   1      -4.195   1.706  13.326  0.50 16.74           C
  """
  h = iotbx.pdb.input(source_info=None, lines=lines).construct_hierarchy()
  asc = h.atom_selection_cache()
  # case 1
  isel = asc.selection("chain A resseq 1 and altloc a").iselection()
  assert isel.size()==1
  h0 = h.select(isel)
  assert approx_equal(h0.atoms().extract_xyz()[0][0], -4.195)
  # case 2
  isel = asc.selection("chain A resseq 1 and altloc A").iselection()
  assert isel.size()==1
  h0 = h.select(isel)
  assert approx_equal(h0.atoms().extract_xyz()[0][0], -4.890)
  ###
  # case 1
  isel = asc.selection("chain A resseq 1 and altid a").iselection()
  assert isel.size()==1
  h0 = h.select(isel)
  assert approx_equal(h0.atoms().extract_xyz()[0][0], -4.195)
  # case 2
  isel = asc.selection("chain A resseq 1 and altid A").iselection()
  assert isel.size()==1
  h0 = h.select(isel)
  assert approx_equal(h0.atoms().extract_xyz()[0][0], -4.890)

if __name__=='__main__':
  test_get_clean_selection_string()
  test_selection_string_from_selection()
  test_selection_string_from_selection2()
  test_avoid_chain_selection()
  test_avoid_chain_selection2()
  test_avoid_hoh()
  test_include_hoh()
  test_selection_with_alternative_conformers()
  test_insertions()
  test_2()
  test_3()
  test_4()
  test_5()
  test_6()
  test_7()
  test_8()
  test_11()
  test_12()
  test_13()
  test_14()
  print("OK")
