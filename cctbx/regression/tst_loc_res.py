from __future__ import absolute_import, division, print_function

from libtbx.test_utils import approx_equal
import iotbx.pdb
from cctbx import maptbx
from  cctbx import maptbx
from six.moves import zip
from iotbx import map_model_manager
from iotbx import map_manager

raw_records = """\
CRYST1   26.880   33.600   31.920  90.00  90.00  90.00 P 1
ATOM      1  N   META1   1       6.457   5.962   7.411  1.00945.14           N
ATOM      2  CA  META1   1       6.950   5.877   5.993  1.00928.61           C
ATOM      3  C   META1   1       8.223   6.699   5.793  1.00910.08           C
ATOM      4  O   META1   1       9.243   6.182   5.328  1.00910.73           O
ATOM      5  CB  META1   1       5.876   6.358   5.004  1.00932.69           C
ATOM      6  CG  META1   1       4.568   5.569   5.015  1.00937.71           C
ATOM      7  SD  META1   1       4.706   3.879   4.402  1.00945.00           S
ATOM      8  CE  META1   1       3.085   3.239   4.820  1.00941.13           C
ATOM      9  N   GLNA1   2       8.145   7.976   6.159  1.00985.67           N
ATOM     10  CA  GLNA1   2       9.224   8.943   5.958  1.00979.66           C
ATOM     11  C   GLNA1   2      10.241   9.025   7.112  1.00974.69           C
ATOM     12  O   GLNA1   2      11.358   9.495   6.896  1.00979.07           O
ATOM     13  CB  GLNA1   2       8.595  10.314   5.711  1.00980.23           C
ATOM     14  CG  GLNA1   2       9.537  11.423   5.259  1.00989.23           C
ATOM     15  CD  GLNA1   2      10.274  11.097   3.976  1.00902.03           C
ATOM     16  OE1 GLNA1   2       9.763  10.385   3.110  1.00904.39           O
ATOM     17  NE2 GLNA1   2      11.484  11.624   3.845  1.00937.18           N
ATOM     18  N   ARGA1   3       9.853   8.577   8.312  1.00964.30           N
ATOM     19  CA  ARGA1   3      10.716   8.560   9.520  1.00966.26           C
ATOM     20  C   ARGA1   3      10.999   9.975  10.053  1.00966.52           C
ATOM     21  O   ARGA1   3      11.841  10.691   9.509  1.00968.91           O
ATOM     22  CB  ARGA1   3      12.026   7.778   9.276  1.00972.03           C
ATOM     23  CG  ARGA1   3      13.003   7.724  10.446  1.00980.32           C
ATOM     24  CD  ARGA1   3      12.400   7.099  11.694  1.00993.10           C
ATOM     25  NE  ARGA1   3      13.347   7.098  12.809  1.00910.27           N
ATOM     26  CZ  ARGA1   3      13.072   6.692  14.050  1.00920.50           C
ATOM     27  NH1 ARGA1   3      14.021   6.744  14.983  1.00930.56           N
ATOM     28  NH2 ARGA1   3      11.860   6.234  14.373  1.00927.11           N
ATOM     29  N   SERA1   4      10.310  10.351  11.133  1.00970.10           N
ATOM     30  CA  SERA1   4      10.354  11.730  11.665  1.00969.72           C
ATOM     31  C   SERA1   4      11.636  12.037  12.458  1.00972.14           C
ATOM     32  O   SERA1   4      12.139  11.164  13.169  1.00972.46           O
ATOM     33  CB  SERA1   4       9.136  12.007  12.554  1.00967.95           C
ATOM     34  OG  SERA1   4       9.282  11.419  13.837  1.00966.18           O
ATOM     35  N   PROA1   5      12.143  13.288  12.370  1.00976.06           N
ATOM     36  CA  PROA1   5      13.432  13.642  12.995  1.00970.63           C
ATOM     37  C   PROA1   5      13.413  13.900  14.515  1.00962.60           C
ATOM     38  O   PROA1   5      14.460  13.795  15.158  1.00953.71           O
ATOM     39  CB  PROA1   5      13.835  14.912  12.244  1.00971.75           C
ATOM     40  CG  PROA1   5      12.538  15.556  11.904  1.00970.02           C
ATOM     41  CD  PROA1   5      11.546  14.450  11.681  1.00969.15           C
ATOM     42  N   VALA1   6      12.247  14.236  15.074  1.00957.91           N
ATOM     43  CA  VALA1   6      12.086  14.442  16.521  1.00956.66           C
ATOM     44  C   VALA1   6      12.466  13.191  17.316  1.00964.35           C
ATOM     45  O   VALA1   6      12.902  13.293  18.463  1.00963.11           O
ATOM     46  CB  VALA1   6      10.646  14.895  16.873  1.00953.86           C
ATOM     47  CG1 VALA1   6       9.629  13.783  16.632  1.00955.04           C
ATOM     48  CG2 VALA1   6      10.564  15.418  18.305  1.00951.21           C
ATOM     49  N   GLUA1   7      12.275  12.023  16.699  1.00981.64           N
ATOM     50  CA  GLUA1   7      12.746  10.740  17.244  1.00986.88           C
ATOM     51  C   GLUA1   7      14.251  10.751  17.542  1.00984.08           C
ATOM     52  O   GLUA1   7      14.682  10.288  18.602  1.00981.69           O
ATOM     53  CB  GLUA1   7      12.428   9.594  16.274  1.00992.25           C
ATOM     54  CG  GLUA1   7      10.942   9.304  16.075  1.00998.55           C
ATOM     55  CD  GLUA1   7      10.255   8.692  17.292  1.00908.01           C
ATOM     56  OE1 GLUA1   7      10.945   8.223  18.225  1.00917.52           O
ATOM     57  OE2 GLUA1   7       9.006   8.674  17.310  1.00911.20           O
ATOM     58  N   ASPA1   8      15.029  11.281  16.597  1.00981.90           N
ATOM     59  CA  ASPA1   8      16.486  11.412  16.728  1.00973.55           C
ATOM     60  C   ASPA1   8      16.859  12.874  16.984  1.00968.82           C
ATOM     61  O   ASPA1   8      17.526  13.518  16.167  1.00967.38           O
ATOM     62  CB  ASPA1   8      17.181  10.885  15.466  1.00972.03           C
ATOM     63  CG  ASPA1   8      16.797   9.452  15.138  1.00970.76           C
ATOM     64  OD1 ASPA1   8      16.670   8.634  16.072  1.00971.90           O
ATOM     65  OD2 ASPA1   8      16.624   9.145  13.939  1.00969.04           O
ATOM     66  N   ALAA1   9      16.415  13.379  18.132  1.00961.64           N
ATOM     67  CA  ALAA1   9      16.603  14.784  18.516  1.00966.44           C
ATOM     68  C   ALAA1   9      17.351  14.937  19.844  1.00967.19           C
ATOM     69  O   ALAA1   9      17.283  14.063  20.712  1.00961.82           O
ATOM     70  CB  ALAA1   9      15.256  15.486  18.591  1.00968.54           C
ATOM     71  N   ASNA1  10      18.056  16.065  19.975  1.00973.18           N
ATOM     72  CA  ASNA1  10      18.795  16.429  21.183  1.00964.26           C
ATOM     73  C   ASNA1  10      18.003  17.403  22.070  1.00968.71           C
ATOM     74  O   ASNA1  10      17.219  18.228  21.577  1.00959.02           O
ATOM     75  CB  ASNA1  10      20.144  17.045  20.791  1.00956.30           C
ATOM     76  CG  ASNA1  10      21.141  17.087  21.939  1.00955.13           C
ATOM     77  OD1 ASNA1  10      20.945  16.476  22.991  1.00949.89           O
ATOM     78  ND2 ASNA1  10      22.228  17.814  21.732  1.00958.95           N
ATOM     79  N   CYSA1  11      18.240  17.302  23.378  1.00 20.00           N
ATOM     80  CA  CYSA1  11      17.509  18.064  24.406  1.00 20.00           C
ATOM     81  C   CYSA1  11      17.457  19.576  24.165  1.00 20.00           C
ATOM     82  O   CYSA1  11      16.456  20.218  24.484  1.00 20.00           O
ATOM     83  CB  CYSA1  11      18.093  17.773  25.796  1.00 20.00           C
ATOM     84  SG  CYSA1  11      19.900  17.783  25.877  1.00 20.00           S
ATOM     85  N   LEUA1  12      18.530  20.132  23.607  1.00 20.00           N
ATOM     86  CA  LEUA1  12      18.578  21.549  23.243  1.00 20.00           C
ATOM     87  C   LEUA1  12      17.781  21.794  21.968  1.00 20.00           C
ATOM     88  O   LEUA1  12      16.997  22.747  21.888  1.00 20.00           O
ATOM     89  CB  LEUA1  12      20.023  22.024  23.056  1.00 20.00           C
ATOM     90  CG  LEUA1  12      20.775  22.412  24.332  1.00 20.00           C
ATOM     91  CD1 LEUA1  12      21.100  21.195  25.188  1.00 20.00           C
ATOM     92  CD2 LEUA1  12      22.051  23.166  23.986  1.00 20.00           C
ATOM     93  N   SERA1  13      17.979  20.922  20.982  1.00 20.00           N
ATOM     94  CA  SERA1  13      17.342  21.072  19.672  1.00 20.00           C
ATOM     95  C   SERA1  13      15.823  20.970  19.753  1.00 20.00           C
ATOM     96  O   SERA1  13      15.111  21.701  19.049  1.00 20.00           O
ATOM     97  CB  SERA1  13      17.897  20.057  18.667  1.00 20.00           C
ATOM     98  OG  SERA1  13      17.826  18.736  19.165  1.00 20.00           O
ATOM     99  N   ARGA1  14      15.340  20.073  20.614  1.00 20.00           N
ATOM    100  CA  ARGA1  14      13.895  19.927  20.840  1.00 20.00           C
ATOM    101  C   ARGA1  14      13.283  21.099  21.632  1.00 20.00           C
ATOM    102  O   ARGA1  14      12.119  21.445  21.419  1.00 20.00           O
ATOM    103  CB  ARGA1  14      13.570  18.579  21.500  1.00 20.00           C
ATOM    104  CG  ARGA1  14      14.072  18.405  22.922  1.00 20.00           C
ATOM    105  CD  ARGA1  14      13.841  16.988  23.418  1.00 20.00           C
ATOM    106  NE  ARGA1  14      12.423  16.714  23.658  1.00 20.00           N
ATOM    107  CZ  ARGA1  14      11.935  15.592  24.195  1.00 20.00           C
ATOM    108  NH1 ARGA1  14      10.619  15.463  24.369  1.00 20.00           N
ATOM    109  NH2 ARGA1  14      12.742  14.594  24.564  1.00 20.00           N
ATOM    110  N   TYRA1  15      14.065  21.700  22.532  1.00 20.00           N
ATOM    111  CA  TYRA1  15      13.624  22.854  23.325  1.00 20.00           C
ATOM    112  C   TYRA1  15      13.416  24.094  22.457  1.00 20.00           C
ATOM    113  O   TYRA1  15      12.343  24.698  22.482  1.00 20.00           O
ATOM    114  CB  TYRA1  15      14.648  23.155  24.422  1.00 20.00           C
ATOM    115  CG  TYRA1  15      14.241  24.242  25.389  1.00 20.00           C
ATOM    116  CD1 TYRA1  15      13.524  23.941  26.547  1.00 20.00           C
ATOM    117  CD2 TYRA1  15      14.590  25.574  25.161  1.00 20.00           C
ATOM    118  CE1 TYRA1  15      13.156  24.935  27.445  1.00 20.00           C
ATOM    119  CE2 TYRA1  15      14.226  26.577  26.050  1.00 20.00           C
ATOM    120  CZ  TYRA1  15      13.510  26.254  27.191  1.00 20.00           C
ATOM    121  OH  TYRA1  15      13.147  27.244  28.077  1.00 20.00           O
ATOM    122  N   PHEA1  16      14.448  24.460  21.697  1.00 20.00           N
ATOM    123  CA  PHEA1  16      14.402  25.621  20.779  1.00 20.00           C
ATOM    124  C   PHEA1  16      13.776  25.358  19.393  1.00 20.00           C
ATOM    125  O   PHEA1  16      13.545  26.313  18.643  1.00 20.00           O
ATOM    126  CB  PHEA1  16      15.814  26.192  20.575  1.00 20.00           C
ATOM    127  CG  PHEA1  16      16.319  26.997  21.734  1.00 20.00           C
ATOM    128  CD1 PHEA1  16      16.018  28.349  21.840  1.00 20.00           C
ATOM    129  CD2 PHEA1  16      17.110  26.410  22.712  1.00 20.00           C
ATOM    130  CE1 PHEA1  16      16.490  29.099  22.909  1.00 20.00           C
ATOM    131  CE2 PHEA1  16      17.581  27.153  23.785  1.00 20.00           C
ATOM    132  CZ  PHEA1  16      17.272  28.501  23.883  1.00 20.00           C
ATOM    133  N   PHEA1  17      13.499  24.089  19.068  1.00 20.00           N
ATOM    134  CA  PHEA1  17      13.061  23.671  17.724  1.00 20.00           C
ATOM    135  C   PHEA1  17      14.112  23.971  16.634  1.00 20.00           C
ATOM    136  O   PHEA1  17      13.961  24.927  15.864  1.00 20.00           O
ATOM    137  CB  PHEA1  17      11.696  24.288  17.357  1.00 20.00           C
ATOM    138  CG  PHEA1  17      10.643  24.113  18.413  1.00 20.00           C
ATOM    139  CD1 PHEA1  17       9.834  22.986  18.423  1.00 20.00           C
ATOM    140  CD2 PHEA1  17      10.451  25.082  19.394  1.00 20.00           C
ATOM    141  CE1 PHEA1  17       8.855  22.825  19.394  1.00 20.00           C
ATOM    142  CE2 PHEA1  17       9.476  24.925  20.369  1.00 20.00           C
ATOM    143  CZ  PHEA1  17       8.675  23.795  20.370  1.00 20.00           C
ATOM    144  N   TRPA1  18      15.177  23.164  16.593  1.00 20.00           N
ATOM    145  CA  TRPA1  18      16.173  23.217  15.500  1.00 20.00           C
ATOM    146  C   TRPA1  18      15.755  22.365  14.287  1.00 20.00           C
ATOM    147  O   TRPA1  18      16.241  22.588  13.175  1.00 20.00           O
ATOM    148  CB  TRPA1  18      17.555  22.750  15.985  1.00 20.00           C
ATOM    149  CG  TRPA1  18      18.326  23.748  16.815  1.00 20.00           C
ATOM    150  CD1 TRPA1  18      17.963  24.262  18.026  1.00 20.00           C
ATOM    151  CD2 TRPA1  18      19.609  24.321  16.505  1.00 20.00           C
ATOM    152  NE1 TRPA1  18      18.931  25.126  18.488  1.00 20.00           N
ATOM    153  CE2 TRPA1  18      19.949  25.185  17.575  1.00 20.00           C
ATOM    154  CE3 TRPA1  18      20.499  24.197  15.427  1.00 20.00           C
ATOM    155  CZ2 TRPA1  18      21.145  25.922  17.598  1.00 20.00           C
ATOM    156  CZ3 TRPA1  18      21.692  24.933  15.451  1.00 20.00           C
ATOM    157  CH2 TRPA1  18      22.000  25.783  16.532  1.00 20.00           C
ATOM    158  N   TRPA1  19      14.858  21.401  14.508  1.00 20.00           N
ATOM    159  CA  TRPA1  19      14.551  20.352  13.512  1.00 20.00           C
ATOM    160  C   TRPA1  19      13.654  20.762  12.341  1.00 20.00           C
ATOM    161  O   TRPA1  19      13.665  20.095  11.303  1.00 20.00           O
ATOM    162  CB  TRPA1  19      13.982  19.092  14.183  1.00 20.00           C
ATOM    163  CG  TRPA1  19      12.647  19.247  14.852  1.00 20.00           C
ATOM    164  CD1 TRPA1  19      12.423  19.644  16.133  1.00 20.00           C
ATOM    165  CD2 TRPA1  19      11.359  18.966  14.289  1.00 20.00           C
ATOM    166  NE1 TRPA1  19      11.079  19.648  16.401  1.00 20.00           N
ATOM    167  CE2 TRPA1  19      10.400  19.233  15.289  1.00 20.00           C
ATOM    168  CE3 TRPA1  19      10.921  18.522  13.036  1.00 20.00           C
ATOM    169  CZ2 TRPA1  19       9.025  19.072  15.080  1.00 20.00           C
ATOM    170  CZ3 TRPA1  19       9.548  18.359  12.825  1.00 20.00           C
ATOM    171  CH2 TRPA1  19       8.619  18.635  13.845  1.00 20.00           C
ATOM    172  N   THRA1  20      12.893  21.842  12.494  1.00 20.00           N
ATOM    173  CA  THRA1  20      12.024  22.334  11.408  1.00 20.00           C
ATOM    174  C   THRA1  20      12.757  22.960  10.196  1.00 20.00           C
ATOM    175  O   THRA1  20      12.143  23.126   9.131  1.00 20.00           O
ATOM    176  CB  THRA1  20      10.981  23.340  11.942  1.00 20.00           C
ATOM    177  OG1 THRA1  20      11.638  24.387  12.665  1.00 20.00           O
ATOM    178  CG2 THRA1  20       9.986  22.643  12.860  1.00 20.00           C
""".splitlines()

raw_records_2 = """\
CRYST1   26.880   33.600   31.920  90.00  90.00  90.00 P 1
ATOM      1  N   META1   1       6.457   5.962   7.411  1.00200.00           N
ATOM      2  CA  META1   1       6.950   5.877   5.993  1.00200.00           C
ATOM      3  C   META1   1       8.223   6.699   5.793  1.00200.00           C
ATOM      4  O   META1   1       9.243   6.182   5.328  1.00200.00           O
ATOM      5  CB  META1   1       5.876   6.358   5.004  1.00200.00           C
ATOM      6  CG  META1   1       4.568   5.569   5.015  1.00200.00           C
ATOM      7  SD  META1   1       4.706   3.879   4.402  1.00200.00           S
ATOM      8  CE  META1   1       3.085   3.239   4.820  1.00200.00           C
ATOM      9  N   GLNA1   2       8.145   7.976   6.159  1.00200.00           N
ATOM     10  CA  GLNA1   2       9.224   8.943   5.958  1.00200.00           C
ATOM     11  C   GLNA1   2      10.241   9.025   7.112  1.00200.00           C
ATOM     12  O   GLNA1   2      11.358   9.495   6.896  1.00200.00           O
ATOM     13  CB  GLNA1   2       8.595  10.314   5.711  1.00200.00           C
ATOM     14  CG  GLNA1   2       9.537  11.423   5.259  1.00200.00           C
ATOM     15  CD  GLNA1   2      10.274  11.097   3.976  1.00200.00           C
ATOM     16  OE1 GLNA1   2       9.763  10.385   3.110  1.00200.00           O
ATOM     17  NE2 GLNA1   2      11.484  11.624   3.845  1.00200.00           N
ATOM     18  N   ARGA1   3       9.853   8.577   8.312  1.00200.00           N
ATOM     19  CA  ARGA1   3      10.716   8.560   9.520  1.00200.00           C
ATOM     20  C   ARGA1   3      10.999   9.975  10.053  1.00200.00           C
ATOM     21  O   ARGA1   3      11.841  10.691   9.509  1.00200.00           O
ATOM     22  CB  ARGA1   3      12.026   7.778   9.276  1.00200.00           C
ATOM     23  CG  ARGA1   3      13.003   7.724  10.446  1.00200.00           C
ATOM     24  CD  ARGA1   3      12.400   7.099  11.694  1.00200.00           C
ATOM     25  NE  ARGA1   3      13.347   7.098  12.809  1.00200.00           N
ATOM     26  CZ  ARGA1   3      13.072   6.692  14.050  1.00200.00           C
ATOM     27  NH1 ARGA1   3      14.021   6.744  14.983  1.00200.00           N
ATOM     28  NH2 ARGA1   3      11.860   6.234  14.373  1.00200.00           N
ATOM     29  N   SERA1   4      10.310  10.351  11.133  1.00200.00           N
ATOM     30  CA  SERA1   4      10.354  11.730  11.665  1.00200.00           C
ATOM     31  C   SERA1   4      11.636  12.037  12.458  1.00200.00           C
ATOM     32  O   SERA1   4      12.139  11.164  13.169  1.00200.00           O
ATOM     33  CB  SERA1   4       9.136  12.007  12.554  1.00200.00           C
ATOM     34  OG  SERA1   4       9.282  11.419  13.837  1.00200.00           O
ATOM     35  N   PROA1   5      12.143  13.288  12.370  1.00200.00           N
ATOM     36  CA  PROA1   5      13.432  13.642  12.995  1.00200.00           C
ATOM     37  C   PROA1   5      13.413  13.900  14.515  1.00200.00           C
ATOM     38  O   PROA1   5      14.460  13.795  15.158  1.00200.00           O
ATOM     39  CB  PROA1   5      13.835  14.912  12.244  1.00200.00           C
ATOM     40  CG  PROA1   5      12.538  15.556  11.904  1.00200.00           C
ATOM     41  CD  PROA1   5      11.546  14.450  11.681  1.00200.00           C
ATOM     42  N   VALA1   6      12.247  14.236  15.074  1.00200.00           N
ATOM     43  CA  VALA1   6      12.086  14.442  16.521  1.00200.00           C
ATOM     44  C   VALA1   6      12.466  13.191  17.316  1.00200.00           C
ATOM     45  O   VALA1   6      12.902  13.293  18.463  1.00200.00           O
ATOM     46  CB  VALA1   6      10.646  14.895  16.873  1.00200.00           C
ATOM     47  CG1 VALA1   6       9.629  13.783  16.632  1.00200.00           C
ATOM     48  CG2 VALA1   6      10.564  15.418  18.305  1.00200.00           C
ATOM     49  N   GLUA1   7      12.275  12.023  16.699  1.00200.00           N
ATOM     50  CA  GLUA1   7      12.746  10.740  17.244  1.00200.00           C
ATOM     51  C   GLUA1   7      14.251  10.751  17.542  1.00200.00           C
ATOM     52  O   GLUA1   7      14.682  10.288  18.602  1.00200.00           O
ATOM     53  CB  GLUA1   7      12.428   9.594  16.274  1.00200.00           C
ATOM     54  CG  GLUA1   7      10.942   9.304  16.075  1.00200.00           C
ATOM     55  CD  GLUA1   7      10.255   8.692  17.292  1.00200.00           C
ATOM     56  OE1 GLUA1   7      10.945   8.223  18.225  1.00200.00           O
ATOM     57  OE2 GLUA1   7       9.006   8.674  17.310  1.00200.00           O
ATOM     58  N   ASPA1   8      15.029  11.281  16.597  1.00200.00           N
ATOM     59  CA  ASPA1   8      16.486  11.412  16.728  1.00200.00           C
ATOM     60  C   ASPA1   8      16.859  12.874  16.984  1.00200.00           C
ATOM     61  O   ASPA1   8      17.526  13.518  16.167  1.00200.00           O
ATOM     62  CB  ASPA1   8      17.181  10.885  15.466  1.00200.00           C
ATOM     63  CG  ASPA1   8      16.797   9.452  15.138  1.00200.00           C
ATOM     64  OD1 ASPA1   8      16.670   8.634  16.072  1.00200.00           O
ATOM     65  OD2 ASPA1   8      16.624   9.145  13.939  1.00200.00           O
ATOM     66  N   ALAA1   9      16.415  13.379  18.132  1.00200.00           N
ATOM     67  CA  ALAA1   9      16.603  14.784  18.516  1.00200.00           C
ATOM     68  C   ALAA1   9      17.351  14.937  19.844  1.00200.00           C
ATOM     69  O   ALAA1   9      17.283  14.063  20.712  1.00200.00           O
ATOM     70  CB  ALAA1   9      15.256  15.486  18.591  1.00200.00           C
ATOM     71  N   ASNA1  10      18.056  16.065  19.975  1.00200.00           N
ATOM     72  CA  ASNA1  10      18.795  16.429  21.183  1.00200.00           C
ATOM     73  C   ASNA1  10      18.003  17.403  22.070  1.00200.00           C
ATOM     74  O   ASNA1  10      17.219  18.228  21.577  1.00200.00           O
ATOM     75  CB  ASNA1  10      20.144  17.045  20.791  1.00200.00           C
ATOM     76  CG  ASNA1  10      21.141  17.087  21.939  1.00200.00           C
ATOM     77  OD1 ASNA1  10      20.945  16.476  22.991  1.00200.00           O
ATOM     78  ND2 ASNA1  10      22.228  17.814  21.732  1.00200.00           N
ATOM     79  N   CYSA1  11      18.240  17.302  23.378  1.00 20.00           N
ATOM     80  CA  CYSA1  11      17.509  18.064  24.406  1.00 20.00           C
ATOM     81  C   CYSA1  11      17.457  19.576  24.165  1.00 20.00           C
ATOM     82  O   CYSA1  11      16.456  20.218  24.484  1.00 20.00           O
ATOM     83  CB  CYSA1  11      18.093  17.773  25.796  1.00 20.00           C
ATOM     84  SG  CYSA1  11      19.900  17.783  25.877  1.00 20.00           S
ATOM     85  N   LEUA1  12      18.530  20.132  23.607  1.00 20.00           N
ATOM     86  CA  LEUA1  12      18.578  21.549  23.243  1.00 20.00           C
ATOM     87  C   LEUA1  12      17.781  21.794  21.968  1.00 20.00           C
ATOM     88  O   LEUA1  12      16.997  22.747  21.888  1.00 20.00           O
ATOM     89  CB  LEUA1  12      20.023  22.024  23.056  1.00 20.00           C
ATOM     90  CG  LEUA1  12      20.775  22.412  24.332  1.00 20.00           C
ATOM     91  CD1 LEUA1  12      21.100  21.195  25.188  1.00 20.00           C
ATOM     92  CD2 LEUA1  12      22.051  23.166  23.986  1.00 20.00           C
ATOM     93  N   SERA1  13      17.979  20.922  20.982  1.00 20.00           N
ATOM     94  CA  SERA1  13      17.342  21.072  19.672  1.00 20.00           C
ATOM     95  C   SERA1  13      15.823  20.970  19.753  1.00 20.00           C
ATOM     96  O   SERA1  13      15.111  21.701  19.049  1.00 20.00           O
ATOM     97  CB  SERA1  13      17.897  20.057  18.667  1.00 20.00           C
ATOM     98  OG  SERA1  13      17.826  18.736  19.165  1.00 20.00           O
ATOM     99  N   ARGA1  14      15.340  20.073  20.614  1.00 20.00           N
ATOM    100  CA  ARGA1  14      13.895  19.927  20.840  1.00 20.00           C
ATOM    101  C   ARGA1  14      13.283  21.099  21.632  1.00 20.00           C
ATOM    102  O   ARGA1  14      12.119  21.445  21.419  1.00 20.00           O
ATOM    103  CB  ARGA1  14      13.570  18.579  21.500  1.00 20.00           C
ATOM    104  CG  ARGA1  14      14.072  18.405  22.922  1.00 20.00           C
ATOM    105  CD  ARGA1  14      13.841  16.988  23.418  1.00 20.00           C
ATOM    106  NE  ARGA1  14      12.423  16.714  23.658  1.00 20.00           N
ATOM    107  CZ  ARGA1  14      11.935  15.592  24.195  1.00 20.00           C
ATOM    108  NH1 ARGA1  14      10.619  15.463  24.369  1.00 20.00           N
ATOM    109  NH2 ARGA1  14      12.742  14.594  24.564  1.00 20.00           N
ATOM    110  N   TYRA1  15      14.065  21.700  22.532  1.00 20.00           N
ATOM    111  CA  TYRA1  15      13.624  22.854  23.325  1.00 20.00           C
ATOM    112  C   TYRA1  15      13.416  24.094  22.457  1.00 20.00           C
ATOM    113  O   TYRA1  15      12.343  24.698  22.482  1.00 20.00           O
ATOM    114  CB  TYRA1  15      14.648  23.155  24.422  1.00 20.00           C
ATOM    115  CG  TYRA1  15      14.241  24.242  25.389  1.00 20.00           C
ATOM    116  CD1 TYRA1  15      13.524  23.941  26.547  1.00 20.00           C
ATOM    117  CD2 TYRA1  15      14.590  25.574  25.161  1.00 20.00           C
ATOM    118  CE1 TYRA1  15      13.156  24.935  27.445  1.00 20.00           C
ATOM    119  CE2 TYRA1  15      14.226  26.577  26.050  1.00 20.00           C
ATOM    120  CZ  TYRA1  15      13.510  26.254  27.191  1.00 20.00           C
ATOM    121  OH  TYRA1  15      13.147  27.244  28.077  1.00 20.00           O
ATOM    122  N   PHEA1  16      14.448  24.460  21.697  1.00 20.00           N
ATOM    123  CA  PHEA1  16      14.402  25.621  20.779  1.00 20.00           C
ATOM    124  C   PHEA1  16      13.776  25.358  19.393  1.00 20.00           C
ATOM    125  O   PHEA1  16      13.545  26.313  18.643  1.00 20.00           O
ATOM    126  CB  PHEA1  16      15.814  26.192  20.575  1.00 20.00           C
ATOM    127  CG  PHEA1  16      16.319  26.997  21.734  1.00 20.00           C
ATOM    128  CD1 PHEA1  16      16.018  28.349  21.840  1.00 20.00           C
ATOM    129  CD2 PHEA1  16      17.110  26.410  22.712  1.00 20.00           C
ATOM    130  CE1 PHEA1  16      16.490  29.099  22.909  1.00 20.00           C
ATOM    131  CE2 PHEA1  16      17.581  27.153  23.785  1.00 20.00           C
ATOM    132  CZ  PHEA1  16      17.272  28.501  23.883  1.00 20.00           C
ATOM    133  N   PHEA1  17      13.499  24.089  19.068  1.00 20.00           N
ATOM    134  CA  PHEA1  17      13.061  23.671  17.724  1.00 20.00           C
ATOM    135  C   PHEA1  17      14.112  23.971  16.634  1.00 20.00           C
ATOM    136  O   PHEA1  17      13.961  24.927  15.864  1.00 20.00           O
ATOM    137  CB  PHEA1  17      11.696  24.288  17.357  1.00 20.00           C
ATOM    138  CG  PHEA1  17      10.643  24.113  18.413  1.00 20.00           C
ATOM    139  CD1 PHEA1  17       9.834  22.986  18.423  1.00 20.00           C
ATOM    140  CD2 PHEA1  17      10.451  25.082  19.394  1.00 20.00           C
ATOM    141  CE1 PHEA1  17       8.855  22.825  19.394  1.00 20.00           C
ATOM    142  CE2 PHEA1  17       9.476  24.925  20.369  1.00 20.00           C
ATOM    143  CZ  PHEA1  17       8.675  23.795  20.370  1.00 20.00           C
ATOM    144  N   TRPA1  18      15.177  23.164  16.593  1.00 20.00           N
ATOM    145  CA  TRPA1  18      16.173  23.217  15.500  1.00 20.00           C
ATOM    146  C   TRPA1  18      15.755  22.365  14.287  1.00 20.00           C
ATOM    147  O   TRPA1  18      16.241  22.588  13.175  1.00 20.00           O
ATOM    148  CB  TRPA1  18      17.555  22.750  15.985  1.00 20.00           C
ATOM    149  CG  TRPA1  18      18.326  23.748  16.815  1.00 20.00           C
ATOM    150  CD1 TRPA1  18      17.963  24.262  18.026  1.00 20.00           C
ATOM    151  CD2 TRPA1  18      19.609  24.321  16.505  1.00 20.00           C
ATOM    152  NE1 TRPA1  18      18.931  25.126  18.488  1.00 20.00           N
ATOM    153  CE2 TRPA1  18      19.949  25.185  17.575  1.00 20.00           C
ATOM    154  CE3 TRPA1  18      20.499  24.197  15.427  1.00 20.00           C
ATOM    155  CZ2 TRPA1  18      21.145  25.922  17.598  1.00 20.00           C
ATOM    156  CZ3 TRPA1  18      21.692  24.933  15.451  1.00 20.00           C
ATOM    157  CH2 TRPA1  18      22.000  25.783  16.532  1.00 20.00           C
ATOM    158  N   TRPA1  19      14.858  21.401  14.508  1.00 20.00           N
ATOM    159  CA  TRPA1  19      14.551  20.352  13.512  1.00 20.00           C
ATOM    160  C   TRPA1  19      13.654  20.762  12.341  1.00 20.00           C
ATOM    161  O   TRPA1  19      13.665  20.095  11.303  1.00 20.00           O
ATOM    162  CB  TRPA1  19      13.982  19.092  14.183  1.00 20.00           C
ATOM    163  CG  TRPA1  19      12.647  19.247  14.852  1.00 20.00           C
ATOM    164  CD1 TRPA1  19      12.423  19.644  16.133  1.00 20.00           C
ATOM    165  CD2 TRPA1  19      11.359  18.966  14.289  1.00 20.00           C
ATOM    166  NE1 TRPA1  19      11.079  19.648  16.401  1.00 20.00           N
ATOM    167  CE2 TRPA1  19      10.400  19.233  15.289  1.00 20.00           C
ATOM    168  CE3 TRPA1  19      10.921  18.522  13.036  1.00 20.00           C
ATOM    169  CZ2 TRPA1  19       9.025  19.072  15.080  1.00 20.00           C
ATOM    170  CZ3 TRPA1  19       9.548  18.359  12.825  1.00 20.00           C
ATOM    171  CH2 TRPA1  19       8.619  18.635  13.845  1.00 20.00           C
ATOM    172  N   THRA1  20      12.893  21.842  12.494  1.00 20.00           N
ATOM    173  CA  THRA1  20      12.024  22.334  11.408  1.00 20.00           C
ATOM    174  C   THRA1  20      12.757  22.960  10.196  1.00 20.00           C
ATOM    175  O   THRA1  20      12.143  23.126   9.131  1.00 20.00           O
ATOM    176  CB  THRA1  20      10.981  23.340  11.942  1.00 20.00           C
ATOM    177  OG1 THRA1  20      11.638  24.387  12.665  1.00 20.00           O
ATOM    178  CG2 THRA1  20       9.986  22.643  12.860  1.00 20.00           C
""".splitlines()


def get_b_and_occ(hierarchy=None,atom_selection=None):
  atom_selection="name CA and (resid 5 or resid 15)"
  asc1=hierarchy.atom_selection_cache()
  sel1 = asc1.selection(string = atom_selection)
  b_values=list(hierarchy.select(sel1).atoms().extract_b())
  resolutions=list(hierarchy.select(sel1).atoms().extract_occ())
  return b_values,resolutions

def make_map_from_pdb(raw_records=None,set_b_iso=None):
    pdb_inp = iotbx.pdb.input(source_info=None, lines=raw_records)
    ph = pdb_inp.construct_hierarchy()
    xrs = pdb_inp.xray_structure_simple()
    if set_b_iso:
      xrs.set_b_iso(value=set_b_iso)
    crystal_symmetry=xrs.crystal_symmetry()
    fc = xrs.structure_factors(d_min = 3, algorithm = "direct").f_calc()
    fft_map = fc.fft_map(resolution_factor = 0.25)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()
    # Create model instead of pdb_hierarchy
    ph.adopt_xray_structure(xrs)
    import mmtbx.model
    model=mmtbx.model.manager(
          model_input=None,
          pdb_hierarchy = ph,
          crystal_symmetry = crystal_symmetry)
    return model,map_data,crystal_symmetry

def tst_0():
  model,map_data,crystal_symmetry=make_map_from_pdb(raw_records=raw_records)
  mm = map_manager.map_manager(
    map_data                   = map_data,
    unit_cell_crystal_symmetry = crystal_symmetry,
    unit_cell_grid             = map_data.all(),
    wrapping                   = False)
  mmm = map_model_manager.map_model_manager(model = model, map_manager = mm)

  method="rscc"
  r  = maptbx.loc_res(
    map_model_manager = mmm,
    method            = method)
  assert  approx_equal([r.result[0], r.result[-1]], [9.7, 3.1], eps=.05)

  method = "fsc"
  r  = maptbx.loc_res(
    map_model_manager = mmm,
    method            = method)
  assert  approx_equal(sorted([round(r,1) for r in set(r.result)]), [1.9, 5.2], eps=.05)

def tst_1():
  for set_b_iso,expected_b_values,expected_occs in zip(
     [-200,0,200],
     [[200.0, 200.0],[-20.0, -20.0],[-160.0, -180.0]],
     [[3.0, 3.0],[3.4, 3.0],[3.2, 3.1]]   ):

    model,map_data,crystal_symmetry=make_map_from_pdb(raw_records=raw_records_2,
     set_b_iso=set_b_iso)

    mm = map_manager.map_manager(
      map_data                   = map_data,
      unit_cell_crystal_symmetry = crystal_symmetry,
      unit_cell_grid             = map_data.all(),
      wrapping                   = False)
    mmm = map_model_manager.map_model_manager(model = model, map_manager = mm)

    method = "rscc_d_min_b"
    r  = maptbx.loc_res(
      map_model_manager = mmm,
      method            = method,
      b_min             = -1000,
      b_max             = 1000,
      b_step            = 1)
    print("result:", list(set(r.result)), "expected:", expected_occs)

if (__name__ == "__main__"):
  tst_0()
  tst_1()
  print("OK")
