from __future__ import absolute_import, division, print_function
from mmtbx.tls import tools
import time
import iotbx.pdb
import mmtbx.tls.tools
from scitbx import matrix
from mmtbx.tls import analysis
import math
from scitbx.array_family import flex
from six.moves import range

pdb_str_CA = """
REMARK   3  TLS DETAILS.
REMARK   3   NUMBER OF TLS GROUPS: 1
REMARK   3   ORIGIN: CENTER OF MASS
REMARK   3   TLS GROUP : 1
REMARK   3    SELECTION: all
REMARK   3    ORIGIN FOR THE GROUP (A):   8.7338  28.3021  16.6793
REMARK   3    T TENSOR
REMARK   3      T11:   0.0695 T22:   0.0842
REMARK   3      T33:   0.0992 T12:   0.0039
REMARK   3      T13:  -0.0041 T23:   0.0014
REMARK   3    L TENSOR
REMARK   3      L11:   2.1725 L22:   0.5881
REMARK   3      L33:   1.1562 L12:   0.4879
REMARK   3      L13:  -0.0881 L23:   0.2404
REMARK   3    S TENSOR
REMARK   3      S11:   0.0512 S12:   0.0543 S13:   0.0992
REMARK   3      S21:  -0.0892 S22:  -0.0659 S23:  -0.0441
REMARK   3      S31:  -0.0042 S32:  -0.0302 S33:   0.0147
CRYST1   35.050   40.500   42.370  90.00  90.00  90.00 P 21 21 21
SCALE1      0.028531  0.000000  0.000000        0.00000
SCALE2      0.000000  0.024691  0.000000        0.00000
SCALE3      0.000000  0.000000  0.023602        0.00000
ATOM      1  CA  THR A   6       6.096  14.546  15.382  1.00  1.21           C
ATOM      2  CA  THR A   7       4.643  17.379  17.472  1.00  0.67           C
ATOM      3  CA  TYR A   8       7.031  20.308  17.626  1.00 -0.75           C
ATOM      4  CA  LYS A   9       6.715  23.191  20.108  1.00 -0.06           C
ATOM      5  CA  LEU A  10       7.484  26.896  19.736  1.00 -0.50           C
ATOM      6  CA  VAL A  11       8.313  28.973  22.814  1.00 -0.68           C
ATOM      7  CA  ILE A  12       7.745  32.666  21.951  1.00  0.00           C
ATOM      8  CA  ASN A  13       9.367  35.454  23.951  1.00  2.73           C
ATOM      9  CA  GLY A  14       8.223  38.436  21.869  1.00  4.26           C
ATOM     10  CA  LYS A  15       7.732  42.071  22.635  1.00  4.20           C
ATOM     11  CA  THR A  16       3.910  41.794  22.427  1.00 -1.15           C
ATOM     12  CA  LEU A  17       3.328  38.080  21.651  1.00 -2.50           C
ATOM     13  CA  LYS A  18       4.268  35.598  24.351  1.00 -2.00           C
ATOM     14  CA  GLY A  19       3.567  31.957  25.117  1.00 -1.89           C
ATOM     15  CA  GLU A  20       3.694  28.633  23.293  1.00 -1.20           C
ATOM     16  CA  THR A  21       2.207  27.042  20.183  1.00 -0.48           C
ATOM     17  CA  THR A  22       2.609  23.707  18.385  1.00 -0.37           C
ATOM     18  CA  THR A  23       2.682  22.183  14.916  1.00  1.08           C
ATOM     19  CA  LYS A  24       2.823  18.625  13.556  1.00  2.05           C
ATOM     20  CA  ALA A  25       5.312  17.719  10.843  1.00 -0.46           C
ATOM     21  CA  VAL A  26       7.260  14.854   9.299  1.00 -1.58           C
ATOM     22  CA  ASP A  27      10.612  16.526  10.001  1.00 -1.52           C
ATOM     23  CA  ALA A  28      12.142  19.541  11.794  1.00  0.14           C
ATOM     24  CA  GLU A  29      12.712  21.542   8.632  1.00  0.11           C
ATOM     25  CA  THR A  30       9.005  21.472   7.817  1.00 -1.33           C
ATOM     26  CA  ALA A  31       8.121  22.548  11.353  1.00  0.02           C
ATOM     27  CA  GLU A  32      10.719  25.317  11.217  1.00  0.13           C
ATOM     28  CA  LYS A  33       9.139  26.736   8.061  1.00  0.27           C
ATOM     29  CA  ALA A  34       5.654  26.682   9.658  1.00  0.18           C
ATOM     30  CA  PHE A  35       6.940  28.471  12.774  1.00 -0.71           C
ATOM     31  CA  LYS A  36       8.915  31.065  10.827  1.00 -0.16           C
ATOM     32  CA  GLN A  37       5.784  31.817   8.802  1.00 -0.29           C
ATOM     33  CA  TYR A  38       3.757  32.125  12.022  1.00 -0.88           C
ATOM     34  CA  ALA A  39       6.322  34.523  13.495  1.00 -0.11           C
ATOM     35  CA  ASN A  40       6.375  36.639  10.306  1.00  0.84           C
ATOM     36  CA  ASP A  41       2.536  36.702  10.145  1.00 -0.52           C
ATOM     37  CA  ASN A  42       2.519  38.093  13.718  1.00 -0.69           C
ATOM     38  CA  GLY A  43       5.296  40.638  13.264  1.00  0.34           C
ATOM     39  CA  VAL A  44       7.921  38.823  15.338  1.00  0.90           C
ATOM     40  CA  ASP A  45      11.526  39.426  14.290  1.00  3.57           C
ATOM     41  CA  GLY A  46      13.865  37.643  16.623  1.00  2.59           C
ATOM     42  CA  VAL A  47      16.614  35.096  17.081  1.00 -1.03           C
ATOM     43  CA  TRP A  48      16.026  31.364  17.129  1.00 -1.89           C
ATOM     44  CA  THR A  49      17.146  28.133  18.719  1.00 -1.92           C
ATOM     45  CA  TYR A  50      16.165  24.543  17.999  1.00 -1.32           C
ATOM     46  CA  ASP A  51      16.716  21.634  20.410  1.00 -0.77           C
ATOM     47  CA  ASP A  52      16.413  18.270  18.691  1.00  0.97           C
ATOM     48  CA  ALA A  53      16.219  16.432  22.049  1.00  4.22           C
ATOM     49  CA  THR A  54      12.971  18.177  23.062  1.00  1.37           C
ATOM     50  CA  LYS A  55      11.609  19.022  19.577  1.00 -0.72           C
ATOM     51  CA  THR A  56      11.365  22.651  20.684  1.00 -0.98           C
ATOM     52  CA  PHE A  57      12.052  25.886  18.846  1.00 -1.21           C
ATOM     53  CA  THR A  58      12.461  29.177  20.755  1.00 -1.19           C
ATOM     54  CA  VAL A  59      12.275  32.720  19.421  1.00 -1.14           C
ATOM     55  CA  THR A  60      13.202  35.851  21.380  1.00 -0.88           C
ATOM     56  CA  GLU A  61      12.989  39.521  20.535  1.00  1.06           C
TER
END
"""

pdb_str_CACON = """
REMARK   3  TLS DETAILS.
REMARK   3   NUMBER OF TLS GROUPS: 1
REMARK   3   ORIGIN: CENTER OF MASS
REMARK   3   TLS GROUP : 1
REMARK   3    SELECTION: all
REMARK   3    ORIGIN FOR THE GROUP (A):   8.6859  28.4733  16.7346
REMARK   3    T TENSOR
REMARK   3      T11:   0.0691 T22:   0.0854
REMARK   3      T33:   0.1047 T12:   0.0048
REMARK   3      T13:  -0.0039 T23:  -0.0005
REMARK   3    L TENSOR
REMARK   3      L11:   2.1972 L22:   0.5110
REMARK   3      L33:   1.3269 L12:   0.5440
REMARK   3      L13:   0.0026 L23:   0.2089
REMARK   3    S TENSOR
REMARK   3      S11:   0.0612 S12:   0.0383 S13:   0.0986
REMARK   3      S21:  -0.0977 S22:  -0.0764 S23:  -0.0418
REMARK   3      S31:  -0.0211 S32:  -0.0272 S33:   0.0152
CRYST1   35.050   40.500   42.370  90.00  90.00  90.00 P 21 21 21
SCALE1      0.028531  0.000000  0.000000        0.00000
SCALE2      0.000000  0.024691  0.000000        0.00000
SCALE3      0.000000  0.000000  0.023602        0.00000
ATOM      1  N   THR A   6       5.628  14.116  14.075  1.00 -0.14           N
ATOM      2  CA  THR A   6       6.096  14.546  15.382  1.00  0.86           C
ATOM      3  C   THR A   6       5.411  15.882  15.727  1.00  0.13           C
ATOM      4  O   THR A   6       5.231  16.716  14.846  1.00  1.04           O
ATOM      5  N   THR A   7       5.084  16.089  17.004  1.00  0.55           N
ATOM      6  CA  THR A   7       4.643  17.379  17.472  1.00  0.42           C
ATOM      7  C   THR A   7       5.853  18.199  17.897  1.00  0.02           C
ATOM      8  O   THR A   7       6.635  17.775  18.761  1.00  3.64           O
ATOM      9  N   TYR A   8       5.998  19.354  17.311  1.00 -1.13           N
ATOM     10  CA  TYR A   8       7.031  20.308  17.626  1.00 -0.96           C
ATOM     11  C   TYR A   8       6.384  21.477  18.386  1.00 -1.13           C
ATOM     12  O   TYR A   8       5.253  21.855  18.136  1.00  1.22           O
ATOM     13  N   LYS A   9       7.151  22.068  19.294  1.00 -0.54           N
ATOM     14  CA  LYS A   9       6.715  23.191  20.108  1.00 -0.24           C
ATOM     15  C   LYS A   9       7.502  24.425  19.761  1.00 -0.91           C
ATOM     16  O   LYS A   9       8.703  24.345  19.462  1.00  0.55           O
ATOM     17  N   LEU A  10       6.854  25.585  19.831  1.00 -0.60           N
ATOM     18  CA  LEU A  10       7.484  26.896  19.736  1.00 -0.67           C
ATOM     19  C   LEU A  10       7.237  27.653  21.022  1.00 -0.77           C
ATOM     20  O   LEU A  10       6.081  27.852  21.422  1.00  0.79           O
ATOM     21  N   VAL A  11       8.311  28.092  21.657  1.00 -0.91           N
ATOM     22  CA  VAL A  11       8.313  28.973  22.814  1.00 -0.87           C
ATOM     23  C   VAL A  11       8.589  30.397  22.266  1.00 -0.99           C
ATOM     24  O   VAL A  11       9.615  30.605  21.632  1.00  0.68           O
ATOM     25  N   ILE A  12       7.658  31.314  22.493  1.00 -0.10           N
ATOM     26  CA  ILE A  12       7.745  32.666  21.951  1.00 -0.24           C
ATOM     27  C   ILE A  12       7.922  33.678  23.059  1.00 -0.20           C
ATOM     28  O   ILE A  12       7.053  33.835  23.928  1.00  0.88           O
ATOM     29  N   ASN A  13       9.057  34.372  23.043  1.00  1.02           N
ATOM     30  CA  ASN A  13       9.367  35.454  23.951  1.00  2.42           C
ATOM     31  C   ASN A  13       9.451  36.755  23.159  1.00  2.65           C
ATOM     32  O   ASN A  13      10.525  37.295  22.962  1.00  3.78           O
ATOM     33  N   GLY A  14       8.318  37.243  22.670  1.00  3.66           N
ATOM     34  CA  GLY A  14       8.223  38.436  21.869  1.00  3.84           C
ATOM     35  C   GLY A  14       7.950  39.692  22.643  1.00  3.14           C
ATOM     36  O   GLY A  14       7.646  39.596  23.871  1.00  6.24           O
ATOM     37  N   LYS A  15       8.054  40.804  21.975  1.00  3.35           N
ATOM     38  CA  LYS A  15       7.732  42.071  22.635  1.00  3.60           C
ATOM     39  C   LYS A  15       6.257  42.190  23.003  1.00 -0.24           C
ATOM     40  O   LYS A  15       5.977  42.698  24.043  1.00  3.00           O
ATOM     41  N   THR A  16       5.349  41.732  22.127  1.00 -1.64           N
ATOM     42  CA  THR A  16       3.910  41.794  22.427  1.00 -1.81           C
ATOM     43  C   THR A  16       3.272  40.435  22.393  1.00 -3.01           C
ATOM     44  O   THR A  16       2.195  40.299  22.987  1.00  0.66           O
ATOM     45  N   LEU A  17       3.854  39.437  21.763  1.00 -2.63           N
ATOM     46  CA  LEU A  17       3.328  38.080  21.651  1.00 -2.98           C
ATOM     47  C   LEU A  17       4.186  37.184  22.514  1.00 -2.79           C
ATOM     48  O   LEU A  17       5.396  37.109  22.337  1.00  0.26           O
ATOM     49  N   LYS A  18       3.557  36.473  23.453  1.00 -2.77           N
ATOM     50  CA  LYS A  18       4.268  35.598  24.351  1.00 -2.40           C
ATOM     51  C   LYS A  18       3.460  34.329  24.578  1.00 -2.72           C
ATOM     52  O   LYS A  18       2.235  34.355  24.675  1.00 -1.00           O
ATOM     53  N   GLY A  19       4.160  33.211  24.736  1.00 -2.52           N
ATOM     54  CA  GLY A  19       3.567  31.957  25.117  1.00 -2.20           C
ATOM     55  C   GLY A  19       4.126  30.802  24.327  1.00 -2.08           C
ATOM     56  O   GLY A  19       5.355  30.755  24.092  1.00 -0.35           O
ATOM     57  N   GLU A  20       3.289  29.842  23.994  1.00 -2.15           N
ATOM     58  CA  GLU A  20       3.694  28.633  23.293  1.00 -1.43           C
ATOM     59  C   GLU A  20       2.621  28.265  22.255  1.00 -1.41           C
ATOM     60  O   GLU A  20       1.437  28.479  22.482  1.00 -0.28           O
ATOM     61  N   THR A  21       3.079  27.640  21.179  1.00 -0.73           N
ATOM     62  CA  THR A  21       2.207  27.042  20.183  1.00 -0.68           C
ATOM     63  C   THR A  21       2.866  25.745  19.700  1.00 -0.64           C
ATOM     64  O   THR A  21       3.984  25.412  20.101  1.00  0.50           O
ATOM     65  N   THR A  22       2.147  25.001  18.859  1.00 -0.43           N
ATOM     66  CA  THR A  22       2.609  23.707  18.385  1.00 -0.54           C
ATOM     67  C   THR A  22       2.283  23.535  16.910  1.00 -0.03           C
ATOM     68  O   THR A  22       1.334  24.152  16.388  1.00  1.60           O
ATOM     69  N   THR A  23       2.999  22.625  16.268  1.00  0.08           N
ATOM     70  CA  THR A  23       2.682  22.183  14.916  1.00  0.94           C
ATOM     71  C   THR A  23       3.082  20.724  14.791  1.00  0.72           C
ATOM     72  O   THR A  23       3.965  20.253  15.484  1.00  3.94           O
ATOM     73  N   LYS A  24       2.489  20.032  13.836  1.00  1.40           N
ATOM     74  CA  LYS A  24       2.823  18.625  13.556  1.00  1.88           C
ATOM     75  C   LYS A  24       3.526  18.598  12.186  1.00  0.44           C
ATOM     76  O   LYS A  24       3.089  19.260  11.239  1.00  1.80           O
ATOM     77  N   ALA A  25       4.612  17.830  12.094  1.00 -0.69           N
ATOM     78  CA  ALA A  25       5.312  17.719  10.843  1.00 -0.66           C
ATOM     79  C   ALA A  25       6.098  16.407  10.770  1.00 -1.76           C
ATOM     80  O   ALA A  25       6.446  15.811  11.785  1.00 -0.29           O
ATOM     81  N   VAL A  26       6.448  16.040   9.546  1.00 -1.94           N
ATOM     82  CA  VAL A  26       7.260  14.854   9.299  1.00 -1.91           C
ATOM     83  C   VAL A  26       8.714  15.045   9.672  1.00 -2.12           C
ATOM     84  O   VAL A  26       9.432  14.072   9.914  1.00  0.90           O
ATOM     85  N   ASP A  27       9.221  16.273   9.700  1.00 -2.07           N
ATOM     86  CA  ASP A  27      10.612  16.526  10.001  1.00 -1.87           C
ATOM     87  C   ASP A  27      10.736  17.942  10.570  1.00 -1.92           C
ATOM     88  O   ASP A  27       9.811  18.770  10.534  1.00 -0.67           O
ATOM     89  N   ALA A  28      11.911  18.229  11.106  1.00 -0.78           N
ATOM     90  CA  ALA A  28      12.142  19.541  11.794  1.00 -0.15           C
ATOM     91  C   ALA A  28      12.174  20.683  10.816  1.00 -0.64           C
ATOM     92  O   ALA A  28      11.745  21.792  11.167  1.00  1.48           O
ATOM     93  N   GLU A  29      12.656  20.492   9.584  1.00 -0.63           N
ATOM     94  CA  GLU A  29      12.712  21.542   8.632  1.00 -0.14           C
ATOM     95  C   GLU A  29      11.307  22.028   8.263  1.00 -1.18           C
ATOM     96  O   GLU A  29      11.082  23.232   8.087  1.00  0.02           O
ATOM     97  N   THR A  30      10.368  21.120   8.131  1.00 -1.46           N
ATOM     98  CA  THR A  30       9.005  21.472   7.817  1.00 -1.49           C
ATOM     99  C   THR A  30       8.367  22.269   8.951  1.00 -1.23           C
ATOM    100  O   THR A  30       7.682  23.269   8.742  1.00 -0.09           O
ATOM    101  N   ALA A  31       8.601  21.802  10.188  1.00 -0.98           N
ATOM    102  CA  ALA A  31       8.121  22.548  11.353  1.00 -0.13           C
ATOM    103  C   ALA A  31       8.748  23.937  11.415  1.00 -0.59           C
ATOM    104  O   ALA A  31       8.051  24.901  11.739  1.00  0.74           O
ATOM    105  N   GLU A  32      10.048  24.036  11.158  1.00 -0.47           N
ATOM    106  CA  GLU A  32      10.719  25.317  11.217  1.00 -0.03           C
ATOM    107  C   GLU A  32      10.069  26.311  10.244  1.00 -0.35           C
ATOM    108  O   GLU A  32       9.834  27.465  10.577  1.00  0.77           O
ATOM    109  N   LYS A  33       9.765  25.875   9.026  1.00 -0.33           N
ATOM    110  CA  LYS A  33       9.139  26.736   8.061  1.00  0.18           C
ATOM    111  C   LYS A  33       7.767  27.197   8.527  1.00 -0.30           C
ATOM    112  O   LYS A  33       7.400  28.364   8.365  1.00  1.03           O
ATOM    113  N   ALA A  34       6.967  26.300   9.133  1.00 -0.05           N
ATOM    114  CA  ALA A  34       5.654  26.682   9.658  1.00  0.11           C
ATOM    115  C   ALA A  34       5.813  27.704  10.762  1.00 -0.79           C
ATOM    116  O   ALA A  34       5.034  28.686  10.833  1.00  0.64           O
ATOM    117  N   PHE A  35       6.745  27.528  11.666  1.00 -0.78           N
ATOM    118  CA  PHE A  35       6.940  28.471  12.774  1.00 -0.83           C
ATOM    119  C   PHE A  35       7.513  29.793  12.317  1.00 -0.91           C
ATOM    120  O   PHE A  35       7.128  30.833  12.863  1.00 -0.22           O
ATOM    121  N   LYS A  36       8.408  29.812  11.346  1.00 -0.44           N
ATOM    122  CA  LYS A  36       8.915  31.065  10.827  1.00 -0.30           C
ATOM    123  C   LYS A  36       7.793  31.837  10.117  1.00 -0.74           C
ATOM    124  O   LYS A  36       7.745  33.048  10.259  1.00  0.39           O
ATOM    125  N   GLN A  37       6.904  31.157   9.421  1.00 -0.76           N
ATOM    126  CA  GLN A  37       5.784  31.817   8.802  1.00 -0.37           C
ATOM    127  C   GLN A  37       4.860  32.399   9.857  1.00 -1.78           C
ATOM    128  O   GLN A  37       4.396  33.535   9.732  1.00 -0.52           O
ATOM    129  N   TYR A  38       4.554  31.635  10.899  1.00 -1.35           N
ATOM    130  CA  TYR A  38       3.757  32.125  12.022  1.00 -1.03           C
ATOM    131  C   TYR A  38       4.366  33.371  12.623  1.00 -2.04           C
ATOM    132  O   TYR A  38       3.658  34.355  12.914  1.00 -1.22           O
ATOM    133  N   ALA A  39       5.680  33.368  12.879  1.00 -1.26           N
ATOM    134  CA  ALA A  39       6.322  34.523  13.495  1.00 -0.33           C
ATOM    135  C   ALA A  39       6.196  35.741  12.534  1.00 -0.77           C
ATOM    136  O   ALA A  39       5.868  36.832  12.978  1.00  0.43           O
ATOM    137  N   ASN A  40       6.473  35.540  11.269  1.00 -0.25           N
ATOM    138  CA  ASN A  40       6.375  36.639  10.306  1.00  0.62           C
ATOM    139  C   ASN A  40       4.948  37.178  10.237  1.00 -0.63           C
ATOM    140  O   ASN A  40       4.744  38.408  10.195  1.00  0.73           O
ATOM    141  N   ASP A  41       3.959  36.296  10.203  1.00 -1.01           N
ATOM    142  CA  ASP A  41       2.536  36.702  10.145  1.00 -0.75           C
ATOM    143  C   ASP A  41       2.135  37.538  11.355  1.00 -0.27           C
ATOM    144  O   ASP A  41       1.210  38.351  11.240  1.00  3.89           O
ATOM    145  N   ASN A  42       2.798  37.366  12.472  1.00 -1.71           N
ATOM    146  CA  ASN A  42       2.519  38.093  13.718  1.00 -1.05           C
ATOM    147  C   ASN A  42       3.500  39.153  14.001  1.00 -0.84           C
ATOM    148  O   ASN A  42       3.496  39.762  15.071  1.00  1.80           O
ATOM    149  N   GLY A  43       4.380  39.512  13.118  1.00 -0.68           N
ATOM    150  CA  GLY A  43       5.296  40.638  13.264  1.00 -0.09           C
ATOM    151  C   GLY A  43       6.444  40.406  14.185  1.00 -0.46           C
ATOM    152  O   GLY A  43       7.010  41.360  14.689  1.00  3.93           O
ATOM    153  N   VAL A  44       6.830  39.148  14.452  1.00 -0.10           N
ATOM    154  CA  VAL A  44       7.921  38.823  15.338  1.00  0.52           C
ATOM    155  C   VAL A  44       9.187  38.697  14.541  1.00  2.55           C
ATOM    156  O   VAL A  44       9.242  37.875  13.641  1.00  8.31           O
ATOM    157  N   ASP A  45      10.215  39.430  14.941  1.00  1.33           N
ATOM    158  CA  ASP A  45      11.526  39.426  14.290  1.00  3.19           C
ATOM    159  C   ASP A  45      12.590  39.242  15.400  1.00  2.39           C
ATOM    160  O   ASP A  45      13.024  40.161  16.071  1.00  5.21           O
ATOM    161  N   GLY A  46      12.977  38.020  15.568  1.00  2.42           N
ATOM    162  CA  GLY A  46      13.865  37.643  16.623  1.00  2.26           C
ATOM    163  C   GLY A  46      14.903  36.645  16.228  1.00 -0.41           C
ATOM    164  O   GLY A  46      15.180  36.396  15.052  1.00  3.36           O
ATOM    165  N   VAL A  47      15.570  36.098  17.262  1.00 -1.18           N
ATOM    166  CA  VAL A  47      16.614  35.096  17.081  1.00 -1.32           C
ATOM    167  C   VAL A  47      16.094  33.770  17.609  1.00 -2.02           C
ATOM    168  O   VAL A  47      15.264  33.714  18.514  1.00 -0.66           O
ATOM    169  N   TRP A  48      16.607  32.691  17.009  1.00 -2.24           N
ATOM    170  CA  TRP A  48      16.026  31.364  17.129  1.00 -2.13           C
ATOM    171  C   TRP A  48      16.996  30.329  17.636  1.00 -2.80           C
ATOM    172  O   TRP A  48      18.172  30.281  17.232  1.00 -2.00           O
ATOM    173  N   THR A  49      16.488  29.379  18.438  1.00 -1.97           N
ATOM    174  CA  THR A  49      17.146  28.133  18.719  1.00 -2.14           C
ATOM    175  C   THR A  49      16.225  26.967  18.352  1.00 -1.90           C
ATOM    176  O   THR A  49      15.003  27.094  18.322  1.00 -0.77           O
ATOM    177  N   TYR A  50      16.851  25.803  18.146  1.00 -1.29           N
ATOM    178  CA  TYR A  50      16.165  24.543  17.999  1.00 -1.57           C
ATOM    179  C   TYR A  50      16.920  23.500  18.836  1.00 -1.97           C
ATOM    180  O   TYR A  50      18.136  23.407  18.737  1.00 -1.13           O
ATOM    181  N   ASP A  51      16.179  22.734  19.627  1.00 -1.26           N
ATOM    182  CA  ASP A  51      16.716  21.634  20.410  1.00 -1.05           C
ATOM    183  C   ASP A  51      16.100  20.316  19.901  1.00 -1.56           C
ATOM    184  O   ASP A  51      14.901  20.089  20.079  1.00 -0.00           O
ATOM    185  N   ASP A  52      16.901  19.511  19.280  1.00 -0.41           N
ATOM    186  CA  ASP A  52      16.413  18.270  18.691  1.00  0.58           C
ATOM    187  C   ASP A  52      16.029  17.280  19.777  1.00  1.69           C
ATOM    188  O   ASP A  52      15.274  16.349  19.487  1.00  6.79           O
ATOM    189  N   ALA A  53      16.572  17.391  20.979  1.00  1.22           N
ATOM    190  CA  ALA A  53      16.219  16.432  22.049  1.00  3.82           C
ATOM    191  C   ALA A  53      14.808  16.646  22.496  1.00  2.74           C
ATOM    192  O   ALA A  53      14.159  15.682  22.904  1.00 10.40           O
ATOM    193  N   THR A  54      14.305  17.841  22.535  1.00  0.72           N
ATOM    194  CA  THR A  54      12.971  18.177  23.062  1.00  1.09           C
ATOM    195  C   THR A  54      11.998  18.581  21.984  1.00 -0.14           C
ATOM    196  O   THR A  54      10.856  18.904  22.252  1.00  2.59           O
ATOM    197  N   LYS A  55      12.459  18.636  20.724  1.00 -0.87           N
ATOM    198  CA  LYS A  55      11.609  19.022  19.577  1.00 -1.00           C
ATOM    199  C   LYS A  55      11.021  20.391  19.797  1.00 -1.06           C
ATOM    200  O   LYS A  55       9.855  20.639  19.464  1.00  1.55           O
ATOM    201  N   THR A  56      11.822  21.313  20.323  1.00 -1.19           N
ATOM    202  CA  THR A  56      11.365  22.651  20.684  1.00 -1.17           C
ATOM    203  C   THR A  56      12.208  23.723  20.010  1.00 -1.28           C
ATOM    204  O   THR A  56      13.443  23.727  20.120  1.00 -0.50           O
ATOM    205  N   PHE A  57      11.506  24.659  19.368  1.00 -1.32           N
ATOM    206  CA  PHE A  57      12.052  25.886  18.846  1.00 -1.39           C
ATOM    207  C   PHE A  57      11.773  27.014  19.836  1.00 -1.48           C
ATOM    208  O   PHE A  57      10.708  27.005  20.471  1.00 -0.08           O
ATOM    209  N   THR A  58      12.672  27.973  19.913  1.00 -1.33           N
ATOM    210  CA  THR A  58      12.461  29.177  20.755  1.00 -1.37           C
ATOM    211  C   THR A  58      12.794  30.401  19.934  1.00 -1.22           C
ATOM    212  O   THR A  58      13.821  30.406  19.257  1.00  0.75           O
ATOM    213  N   VAL A  59      11.971  31.424  20.012  1.00 -1.18           N
ATOM    214  CA  VAL A  59      12.275  32.720  19.421  1.00 -1.36           C
ATOM    215  C   VAL A  59      12.184  33.794  20.482  1.00 -1.53           C
ATOM    216  O   VAL A  59      11.223  33.845  21.246  1.00  0.21           O
ATOM    217  N   THR A  60      13.188  34.674  20.499  1.00 -1.25           N
ATOM    218  CA  THR A  60      13.202  35.851  21.380  1.00 -1.17           C
ATOM    219  C   THR A  60      13.440  37.095  20.606  1.00 -1.70           C
ATOM    220  O   THR A  60      14.327  37.147  19.737  1.00  0.23           O
ATOM    221  N   GLU A  61      12.727  38.144  20.962  1.00 -0.48           N
ATOM    222  CA  GLU A  61      12.989  39.521  20.535  1.00  0.66           C
ATOM    223  C   GLU A  61      13.742  40.237  21.642  1.00  2.54           C
ATOM    224  O   GLU A  61      14.145  39.578  22.584  1.00  8.39           O
TER
END
"""

def exercise_00(pdb_str, formula):
  """
  CA and main-chain fragments of 2igd.
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  asc = pdb_hierarchy.atom_selection_cache()
  cs = pdb_inp.crystal_symmetry_from_cryst1()
  tls_extract = mmtbx.tls.tools.tls_from_pdb_inp(
    remark_3_records = pdb_inp.extract_remark_iii_records(3),
    pdb_hierarchy    = pdb_hierarchy)
  #
  deg_to_rad_scale = math.pi/180
  tls_params_one_group = tls_extract.tls_params[0]
  T = matrix.sym(sym_mat3=tls_params_one_group.t)
  L = matrix.sym(sym_mat3=tls_params_one_group.l)
  S = matrix.sqr(tls_params_one_group.s)
  origin = tls_params_one_group.origin
  tlso = tools.tlso(
    t      = T.as_sym_mat3(),
    l      = L.as_sym_mat3(),
    s      = S,
    origin = origin)
  log = open("analysis.log","w")
  r = analysis.run(T=T, L=L*(deg_to_rad_scale**2), S=S*deg_to_rad_scale,
    log=log, find_t_S_using_formula=formula).self_check(show=False)
  log.close()
  #
  rs = flex.double()
  for trial in range(10):
    o = tools.u_tls_vs_u_ens(pdb_str=pdb_str,
      dx       = r.dx,
      dy       = r.dy,
      dz       = r.dz,
      sx       = r.sx,
      sy       = r.sy,
      sz       = r.sz,
      lx       = r.l_x,
      ly       = r.l_y,
      lz       = r.l_z,
      tx       = r.tx,
      ty       = r.ty,
      tz       = r.tz,
      vx       = r.v_x,
      vy       = r.v_y,
      vz       = r.v_z,
      w_M_lx   = r.w_M_lx,
      w_M_ly   = r.w_M_ly,
      w_M_lz   = r.w_M_lz,
      origin   = origin,
      n_models = 10000,
      assert_similarity=False)
    rs.append(o.r)
  return flex.mean(rs)

if (__name__ == "__main__"):
  t0 = time.time()
  for formula in ["10","11"]:
    print("formula:", formula)
    for i, pdb_str in enumerate([pdb_str_CA, pdb_str_CACON]):
      r = exercise_00(pdb_str=pdb_str, formula=formula)
      print("  ", i, r)
      if(formula=="10"):
        if(i==0): assert r<0.04
        if(i==1): assert r>0.08
      if(formula=="11"):
        if(i==0): assert r<0.02
        if(i==1): assert r<0.02
