from __future__ import absolute_import, division, print_function
import iotbx.pdb
import time
from cctbx import maptbx
from libtbx.test_utils import approx_equal
import mmtbx.masks

pdb_str_1 = """
CRYST1   13.475   30.191   11.490  90.00  90.00  90.00 P 1
ATOM      1  N   ALA E   1       6.622  27.477   5.274  1.00 50.00           N
ATOM      2  CA  ALA E   1       6.323  27.055   6.635  1.00 50.00           C
ATOM      3  C   ALA E   1       5.442  25.808   6.656  1.00 50.00           C
ATOM      4  O   ALA E   1       5.848  24.796   7.227  1.00 50.00           O
ATOM      5  CB  ALA E   1       5.675  28.191   7.414  1.00 50.00           C
ATOM      6  N   ALA E   2       4.256  25.870   6.012  1.00 50.00           N
ATOM      7  CA  ALA E   2       3.288  24.768   5.916  1.00 50.00           C
ATOM      8  C   ALA E   2       3.847  23.573   5.125  1.00 50.00           C
ATOM      9  O   ALA E   2       3.553  22.422   5.462  1.00 50.00           O
ATOM     10  CB  ALA E   2       2.000  25.267   5.276  1.00 50.00           C
ATOM     11  N   ALA E   3       4.664  23.856   4.084  1.00 50.00           N
ATOM     12  CA  ALA E   3       5.313  22.854   3.231  1.00 50.00           C
ATOM     13  C   ALA E   3       6.391  22.082   3.994  1.00 50.00           C
ATOM     14  O   ALA E   3       6.616  20.906   3.687  1.00 50.00           O
ATOM     15  CB  ALA E   3       5.915  23.517   2.000  1.00 50.00           C
ATOM     16  N   ALA E   4       7.063  22.744   4.976  1.00 50.00           N
ATOM     17  CA  ALA E   4       8.089  22.116   5.822  1.00 50.00           C
ATOM     18  C   ALA E   4       7.441  21.043   6.703  1.00 50.00           C
ATOM     19  O   ALA E   4       7.888  19.897   6.671  1.00 50.00           O
ATOM     20  CB  ALA E   4       8.785  23.160   6.683  1.00 50.00           C
ATOM     21  N   ALA E   5       6.338  21.398   7.423  1.00 50.00           N
ATOM     22  CA  ALA E   5       5.573  20.488   8.298  1.00 50.00           C
ATOM     23  C   ALA E   5       5.042  19.288   7.513  1.00 50.00           C
ATOM     24  O   ALA E   5       4.981  18.179   8.048  1.00 50.00           O
ATOM     25  CB  ALA E   5       4.433  21.233   8.975  1.00 50.00           C
ATOM     26  N   ALA E   6       4.715  19.518   6.221  1.00 50.00           N
ATOM     27  CA  ALA E   6       4.258  18.514   5.260  1.00 50.00           C
ATOM     28  C   ALA E   6       5.440  17.620   4.865  1.00 50.00           C
ATOM     29  O   ALA E   6       5.361  16.405   5.060  1.00 50.00           O
ATOM     30  CB  ALA E   6       3.664  19.193   4.032  1.00 50.00           C
ATOM     31  N   ALA E   7       6.556  18.234   4.374  1.00 50.00           N
ATOM     32  CA  ALA E   7       7.788  17.543   3.962  1.00 50.00           C
ATOM     33  C   ALA E   7       8.411  16.694   5.075  1.00 50.00           C
ATOM     34  O   ALA E   7       8.954  15.625   4.782  1.00 50.00           O
ATOM     35  CB  ALA E   7       8.806  18.540   3.433  1.00 50.00           C
ATOM     36  N   ALA E   8       8.327  17.168   6.343  1.00 50.00           N
ATOM     37  CA  ALA E   8       8.846  16.459   7.515  1.00 50.00           C
ATOM     38  C   ALA E   8       8.009  15.197   7.768  1.00 50.00           C
ATOM     39  O   ALA E   8       8.568  14.097   7.858  1.00 50.00           O
ATOM     40  CB  ALA E   8       8.833  17.368   8.742  1.00 50.00           C
ATOM     41  N   ALA E   9       6.664  15.357   7.803  1.00 50.00           N
ATOM     42  CA  ALA E   9       5.716  14.261   7.994  1.00 50.00           C
ATOM     43  C   ALA E   9       5.754  13.279   6.822  1.00 50.00           C
ATOM     44  O   ALA E   9       5.654  12.075   7.051  1.00 50.00           O
ATOM     45  CB  ALA E   9       4.314  14.805   8.193  1.00 50.00           C
ATOM     46  N   ALA E  10       5.942  13.791   5.578  1.00 50.00           N
ATOM     47  CA  ALA E  10       6.063  12.985   4.357  1.00 50.00           C
ATOM     48  C   ALA E  10       7.322  12.116   4.421  1.00 50.00           C
ATOM     49  O   ALA E  10       7.235  10.910   4.169  1.00 50.00           O
ATOM     50  CB  ALA E  10       6.105  13.881   3.130  1.00 50.00           C
ATOM     51  N   ALA E  11       8.480  12.722   4.803  1.00 50.00           N
ATOM     52  CA  ALA E  11       9.754  12.018   4.967  1.00 50.00           C
ATOM     53  C   ALA E  11       9.631  10.976   6.089  1.00 50.00           C
ATOM     54  O   ALA E  11      10.074   9.837   5.908  1.00 50.00           O
ATOM     55  CB  ALA E  11      10.870  13.004   5.279  1.00 50.00           C
ATOM     56  N   ALA E  12       8.972  11.356   7.217  1.00 50.00           N
ATOM     57  CA  ALA E  12       8.721  10.483   8.371  1.00 50.00           C
ATOM     58  C   ALA E  12       7.851   9.277   7.972  1.00 50.00           C
ATOM     59  O   ALA E  12       8.201   8.140   8.310  1.00 50.00           O
ATOM     60  CB  ALA E  12       8.052  11.269   9.490  1.00 50.00           C
ATOM     61  N   ALA E  13       6.748   9.535   7.215  1.00 50.00           N
ATOM     62  CA  ALA E  13       5.821   8.516   6.708  1.00 50.00           C
ATOM     63  C   ALA E  13       6.530   7.563   5.752  1.00 50.00           C
ATOM     64  O   ALA E  13       6.266   6.360   5.793  1.00 50.00           O
ATOM     65  CB  ALA E  13       4.636   9.173   6.008  1.00 50.00           C
ATOM     66  N   ALA E  14       7.437   8.112   4.903  1.00 50.00           N
ATOM     67  CA  ALA E  14       8.246   7.367   3.939  1.00 50.00           C
ATOM     68  C   ALA E  14       9.277   6.511   4.676  1.00 50.00           C
ATOM     69  O   ALA E  14       9.512   5.363   4.286  1.00 50.00           O
ATOM     70  CB  ALA E  14       8.942   8.329   2.986  1.00 50.00           C
ATOM     71  N   ALA E  15       9.863   7.069   5.767  1.00 50.00           N
ATOM     72  CA  ALA E  15      10.847   6.407   6.629  1.00 50.00           C
ATOM     73  C   ALA E  15      10.213   5.257   7.427  1.00 50.00           C
ATOM     74  O   ALA E  15      10.922   4.322   7.816  1.00 50.00           O
ATOM     75  CB  ALA E  15      11.475   7.418   7.573  1.00 50.00           C
ATOM     76  N   ALA E  16       8.879   5.325   7.652  1.00 50.00           N
ATOM     77  CA  ALA E  16       8.099   4.317   8.370  1.00 50.00           C
ATOM     78  C   ALA E  16       7.822   3.106   7.488  1.00 50.00           C
ATOM     79  O   ALA E  16       7.658   2.000   7.999  1.00 50.00           O
ATOM     80  CB  ALA E  16       6.790   4.918   8.852  1.00 50.00           C
TER
END
"""

pdb_str_2 = """
CRYST1   13.475   30.191   11.490  90.00  90.00  90.00 P 1
ATOM      1  N   ALA E   1      11.622  32.477  10.274  1.00 50.00           N
ATOM      2  CA  ALA E   1      11.323  32.055  11.635  1.00 50.00           C
ATOM      3  C   ALA E   1      10.442  30.808  11.656  1.00 50.00           C
ATOM      4  O   ALA E   1      10.848  29.796  12.227  1.00 50.00           O
ATOM      5  CB  ALA E   1      10.675  33.191  12.414  1.00 50.00           C
ATOM      6  N   ALA E   2       9.256  30.870  11.012  1.00 50.00           N
ATOM      7  CA  ALA E   2       8.288  29.768  10.916  1.00 50.00           C
ATOM      8  C   ALA E   2       8.847  28.573  10.125  1.00 50.00           C
ATOM      9  O   ALA E   2       8.553  27.422  10.462  1.00 50.00           O
ATOM     10  CB  ALA E   2       7.000  30.267  10.276  1.00 50.00           C
ATOM     11  N   ALA E   3       9.664  28.856   9.084  1.00 50.00           N
ATOM     12  CA  ALA E   3      10.313  27.854   8.231  1.00 50.00           C
ATOM     13  C   ALA E   3      11.391  27.082   8.994  1.00 50.00           C
ATOM     14  O   ALA E   3      11.616  25.906   8.687  1.00 50.00           O
ATOM     15  CB  ALA E   3      10.915  28.517   7.000  1.00 50.00           C
ATOM     16  N   ALA E   4      12.063  27.744   9.976  1.00 50.00           N
ATOM     17  CA  ALA E   4      13.089  27.116  10.822  1.00 50.00           C
ATOM     18  C   ALA E   4      12.441  26.043  11.703  1.00 50.00           C
ATOM     19  O   ALA E   4      12.888  24.897  11.671  1.00 50.00           O
ATOM     20  CB  ALA E   4      13.785  28.160  11.683  1.00 50.00           C
ATOM     21  N   ALA E   5      11.338  26.398  12.423  1.00 50.00           N
ATOM     22  CA  ALA E   5      10.573  25.488  13.298  1.00 50.00           C
ATOM     23  C   ALA E   5      10.042  24.288  12.513  1.00 50.00           C
ATOM     24  O   ALA E   5       9.981  23.179  13.048  1.00 50.00           O
ATOM     25  CB  ALA E   5       9.433  26.233  13.975  1.00 50.00           C
ATOM     26  N   ALA E   6       9.715  24.518  11.221  1.00 50.00           N
ATOM     27  CA  ALA E   6       9.258  23.514  10.260  1.00 50.00           C
ATOM     28  C   ALA E   6      10.440  22.620   9.865  1.00 50.00           C
ATOM     29  O   ALA E   6      10.361  21.405  10.060  1.00 50.00           O
ATOM     30  CB  ALA E   6       8.664  24.193   9.032  1.00 50.00           C
ATOM     31  N   ALA E   7      11.556  23.234   9.374  1.00 50.00           N
ATOM     32  CA  ALA E   7      12.788  22.543   8.962  1.00 50.00           C
ATOM     33  C   ALA E   7      13.411  21.694  10.075  1.00 50.00           C
ATOM     34  O   ALA E   7      13.954  20.625   9.782  1.00 50.00           O
ATOM     35  CB  ALA E   7      13.806  23.540   8.433  1.00 50.00           C
ATOM     36  N   ALA E   8      13.327  22.168  11.343  1.00 50.00           N
ATOM     37  CA  ALA E   8      13.846  21.459  12.515  1.00 50.00           C
ATOM     38  C   ALA E   8      13.009  20.197  12.768  1.00 50.00           C
ATOM     39  O   ALA E   8      13.568  19.097  12.858  1.00 50.00           O
ATOM     40  CB  ALA E   8      13.833  22.368  13.742  1.00 50.00           C
ATOM     41  N   ALA E   9      11.664  20.357  12.803  1.00 50.00           N
ATOM     42  CA  ALA E   9      10.716  19.261  12.994  1.00 50.00           C
ATOM     43  C   ALA E   9      10.754  18.279  11.822  1.00 50.00           C
ATOM     44  O   ALA E   9      10.654  17.075  12.051  1.00 50.00           O
ATOM     45  CB  ALA E   9       9.314  19.805  13.193  1.00 50.00           C
ATOM     46  N   ALA E  10      10.942  18.791  10.578  1.00 50.00           N
ATOM     47  CA  ALA E  10      11.063  17.985   9.357  1.00 50.00           C
ATOM     48  C   ALA E  10      12.322  17.116   9.421  1.00 50.00           C
ATOM     49  O   ALA E  10      12.235  15.910   9.169  1.00 50.00           O
ATOM     50  CB  ALA E  10      11.105  18.881   8.130  1.00 50.00           C
ATOM     51  N   ALA E  11      13.480  17.722   9.803  1.00 50.00           N
ATOM     52  CA  ALA E  11      14.754  17.018   9.967  1.00 50.00           C
ATOM     53  C   ALA E  11      14.631  15.976  11.089  1.00 50.00           C
ATOM     54  O   ALA E  11      15.074  14.837  10.908  1.00 50.00           O
ATOM     55  CB  ALA E  11      15.870  18.004  10.279  1.00 50.00           C
ATOM     56  N   ALA E  12      13.972  16.356  12.217  1.00 50.00           N
ATOM     57  CA  ALA E  12      13.721  15.483  13.371  1.00 50.00           C
ATOM     58  C   ALA E  12      12.851  14.277  12.972  1.00 50.00           C
ATOM     59  O   ALA E  12      13.201  13.140  13.310  1.00 50.00           O
ATOM     60  CB  ALA E  12      13.052  16.269  14.490  1.00 50.00           C
ATOM     61  N   ALA E  13      11.748  14.535  12.215  1.00 50.00           N
ATOM     62  CA  ALA E  13      10.821  13.516  11.708  1.00 50.00           C
ATOM     63  C   ALA E  13      11.530  12.563  10.752  1.00 50.00           C
ATOM     64  O   ALA E  13      11.266  11.360  10.793  1.00 50.00           O
ATOM     65  CB  ALA E  13       9.636  14.173  11.008  1.00 50.00           C
ATOM     66  N   ALA E  14      12.437  13.112   9.903  1.00 50.00           N
ATOM     67  CA  ALA E  14      13.246  12.367   8.939  1.00 50.00           C
ATOM     68  C   ALA E  14      14.277  11.511   9.676  1.00 50.00           C
ATOM     69  O   ALA E  14      14.512  10.363   9.286  1.00 50.00           O
ATOM     70  CB  ALA E  14      13.942  13.329   7.986  1.00 50.00           C
ATOM     71  N   ALA E  15      14.863  12.069  10.767  1.00 50.00           N
ATOM     72  CA  ALA E  15      15.847  11.407  11.629  1.00 50.00           C
ATOM     73  C   ALA E  15      15.213  10.257  12.427  1.00 50.00           C
ATOM     74  O   ALA E  15      15.922   9.322  12.816  1.00 50.00           O
ATOM     75  CB  ALA E  15      16.475  12.418  12.573  1.00 50.00           C
ATOM     76  N   ALA E  16      13.879  10.325  12.652  1.00 50.00           N
ATOM     77  CA  ALA E  16      13.099   9.317  13.370  1.00 50.00           C
ATOM     78  C   ALA E  16      12.822   8.106  12.488  1.00 50.00           C
ATOM     79  O   ALA E  16      12.658   7.000  12.999  1.00 50.00           O
ATOM     80  CB  ALA E  16      11.790   9.918  13.852  1.00 50.00           C
TER
END
"""

pdb_str_3 = """
CRYST1   23.475   30.191   21.490  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   ALA E   1      11.658  32.727  10.658  1.00 50.00           N
ATOM      2  CA  ALA E   1      11.151  32.395  11.984  1.00 50.00           C
ATOM      3  C   ALA E   1      10.253  31.164  11.939  1.00 50.00           C
ATOM      4  O   ALA E   1      10.424  30.229  12.721  1.00 50.00           O
ATOM      5  CB  ALA E   1      10.402  33.579  12.578  1.00 50.00           C
ATOM      6  N   ALA E   2       9.294  31.170  11.018  1.00 50.00           N
ATOM      7  CA  ALA E   2       8.368  30.054  10.869  1.00 50.00           C
ATOM      8  C   ALA E   2       8.951  28.970   9.969  1.00 50.00           C
ATOM      9  O   ALA E   2       8.528  27.815  10.019  1.00 50.00           O
ATOM     10  CB  ALA E   2       7.035  30.542  10.321  1.00 50.00           C
ATOM     11  N   ALA E   3       9.923  29.350   9.147  1.00 50.00           N
ATOM     12  CA  ALA E   3      10.566  28.411   8.235  1.00 50.00           C
ATOM     13  C   ALA E   3      11.651  27.608   8.944  1.00 50.00           C
ATOM     14  O   ALA E   3      11.988  26.499   8.529  1.00 50.00           O
ATOM     15  CB  ALA E   3      11.145  29.148   7.037  1.00 50.00           C
ATOM     16  N   ALA E   4      12.196  28.176  10.015  1.00 50.00           N
ATOM     17  CA  ALA E   4      13.244  27.515  10.783  1.00 50.00           C
ATOM     18  C   ALA E   4      12.660  26.464  11.721  1.00 50.00           C
ATOM     19  O   ALA E   4      13.315  25.472  12.042  1.00 50.00           O
ATOM     20  CB  ALA E   4      14.053  28.538  11.566  1.00 50.00           C
ATOM     21  N   ALA E   5      11.425  26.688  12.157  1.00 50.00           N
ATOM     22  CA  ALA E   5      10.751  25.761  13.059  1.00 50.00           C
ATOM     23  C   ALA E   5      10.177  24.570  12.298  1.00 50.00           C
ATOM     24  O   ALA E   5      10.033  23.480  12.851  1.00 50.00           O
ATOM     25  CB  ALA E   5       9.654  26.477  13.833  1.00 50.00           C
ATOM     26  N   ALA E   6       9.852  24.787  11.028  1.00 50.00           N
ATOM     27  CA  ALA E   6       9.294  23.733  10.189  1.00 50.00           C
ATOM     28  C   ALA E   6      10.390  22.825   9.643  1.00 50.00           C
ATOM     29  O   ALA E   6      10.148  21.657   9.338  1.00 50.00           O
ATOM     30  CB  ALA E   6       8.484  24.334   9.051  1.00 50.00           C
ATOM     31  N   ALA E   7      11.597  23.369   9.521  1.00 50.00           N
ATOM     32  CA  ALA E   7      12.732  22.609   9.011  1.00 50.00           C
ATOM     33  C   ALA E   7      13.347  21.738  10.101  1.00 50.00           C
ATOM     34  O   ALA E   7      13.954  20.706   9.815  1.00 50.00           O
ATOM     35  CB  ALA E   7      13.777  23.547   8.425  1.00 50.00           C
ATOM     36  N   ALA E   8      13.185  22.160  11.351  1.00 50.00           N
ATOM     37  CA  ALA E   8      13.723  21.420  12.485  1.00 50.00           C
ATOM     38  C   ALA E   8      12.833  20.234  12.842  1.00 50.00           C
ATOM     39  O   ALA E   8      13.307  19.225  13.364  1.00 50.00           O
ATOM     40  CB  ALA E   8      13.894  22.339  13.685  1.00 50.00           C
ATOM     41  N   ALA E   9      11.542  20.363  12.558  1.00 50.00           N
ATOM     42  CA  ALA E   9      10.583  19.304  12.849  1.00 50.00           C
ATOM     43  C   ALA E   9      10.585  18.243  11.753  1.00 50.00           C
ATOM     44  O   ALA E   9      10.271  17.079  12.001  1.00 50.00           O
ATOM     45  CB  ALA E   9       9.189  19.884  13.027  1.00 50.00           C
ATOM     46  N   ALA E  10      10.940  18.654  10.540  1.00 50.00           N
ATOM     47  CA  ALA E  10      10.983  17.741   9.404  1.00 50.00           C
ATOM     48  C   ALA E  10      12.259  16.906   9.416  1.00 50.00           C
ATOM     49  O   ALA E  10      12.282  15.782   8.914  1.00 50.00           O
ATOM     50  CB  ALA E  10      10.864  18.512   8.099  1.00 50.00           C
ATOM     51  N   ALA E  11      13.320  17.463   9.992  1.00 50.00           N
ATOM     52  CA  ALA E  11      14.601  16.771  10.070  1.00 50.00           C
ATOM     53  C   ALA E  11      14.610  15.753  11.206  1.00 50.00           C
ATOM     54  O   ALA E  11      15.320  14.750  11.147  1.00 50.00           O
ATOM     55  CB  ALA E  11      15.734  17.771  10.241  1.00 50.00           C
ATOM     56  N   ALA E  12      13.817  16.020  12.238  1.00 50.00           N
ATOM     57  CA  ALA E  12      13.732  15.128  13.389  1.00 50.00           C
ATOM     58  C   ALA E  12      12.813  13.945  13.104  1.00 50.00           C
ATOM     59  O   ALA E  12      12.964  12.874  13.690  1.00 50.00           O
ATOM     60  CB  ALA E  12      13.254  15.890  14.615  1.00 50.00           C
ATOM     61  N   ALA E  13      11.860  14.147  12.200  1.00 50.00           N
ATOM     62  CA  ALA E  13      10.915  13.099  11.835  1.00 50.00           C
ATOM     63  C   ALA E  13      11.506  12.164  10.785  1.00 50.00           C
ATOM     64  O   ALA E  13      11.087  11.013  10.658  1.00 50.00           O
ATOM     65  CB  ALA E  13       9.614  13.708  11.334  1.00 50.00           C
ATOM     66  N   ALA E  14      12.482  12.666  10.036  1.00 50.00           N
ATOM     67  CA  ALA E  14      13.132  11.877   8.996  1.00 50.00           C
ATOM     68  C   ALA E  14      14.227  10.991   9.580  1.00 50.00           C
ATOM     69  O   ALA E  14      14.585   9.966   9.000  1.00 50.00           O
ATOM     70  CB  ALA E  14      13.702  12.787   7.919  1.00 50.00           C
ATOM     71  N   ALA E  15      14.755  11.392  10.732  1.00 50.00           N
ATOM     72  CA  ALA E  15      15.809  10.636  11.397  1.00 50.00           C
ATOM     73  C   ALA E  15      15.227   9.611  12.365  1.00 50.00           C
ATOM     74  O   ALA E  15      15.925   8.704  12.818  1.00 50.00           O
ATOM     75  CB  ALA E  15      16.759  11.576  12.124  1.00 50.00           C
ATOM     76  N   ALA E  16      13.944   9.762  12.677  1.00 50.00           N
ATOM     77  CA  ALA E  16      13.265   8.851  13.591  1.00 50.00           C
ATOM     78  C   ALA E  16      12.707   7.641  12.849  1.00 50.00           C
ATOM     79  O   ALA E  16      12.597   6.552  13.413  1.00 50.00           O
ATOM     80  CB  ALA E  16      12.157   9.577  14.337  1.00 50.00           C
TER
"""

def exercise():
  """ Compare simple mask (maptbx/mask.h) with classic mask. Does not account
      for site multiplicity due to symmetry. Also, this is an alternative
      exercise for asu mask.
  """
  for i_pdb, pdb_str in enumerate([pdb_str_1, pdb_str_2, pdb_str_3]):
    for solvent_radius in [0, 1.0]:
      print("file %d solvent_raius: %3.1f"%(i_pdb, solvent_radius))
      xrs = iotbx.pdb.input(source_info = None,
        lines = pdb_str).xray_structure_simple()
      mp = mmtbx.masks.mask_master_params.extract()
      mp.solvent_radius=solvent_radius
      mp.shrink_truncation_radius=0.0
      mp.grid_step_factor=4
      mmtbx_masks_asu_mask_obj = mmtbx.masks.asu_mask(
        xray_structure = xrs.expand_to_p1(sites_mod_positive=True), # Must be P1
        d_min          = 1,
        mask_params    = mp)
      mask_data = mmtbx_masks_asu_mask_obj.mask_data_whole_uc()
      c1o,c0o,s,csf = mask_data.count(1), mask_data.count(0), mask_data.size(),\
        mmtbx_masks_asu_mask_obj.asu_mask.contact_surface_fraction
      assert c1o+c0o==s
      f1o, f0o = c1o*100./s, c0o*100./s
      print("  old", f1o, f0o, csf)
      mask_data = maptbx.mask(
        xray_structure = xrs, # Expanded to P1 internally
        n_real         = mask_data.focus(),
        solvent_radius = solvent_radius)
      c1n, c0n, s = mask_data.count(1), mask_data.count(0), mask_data.size()
      assert c1n+c0n==s
      f1n, f0n = c1n*100./s, c0n*100./s
      print("  new:", f1n, f0n)
      assert approx_equal(f1o, f1n, 1)
      assert approx_equal(f0o, f0n, 1)
      assert approx_equal(f1o, csf*100, 1)

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("Total time: %-8.4f"%(time.time()-t0))
  print("OK")
