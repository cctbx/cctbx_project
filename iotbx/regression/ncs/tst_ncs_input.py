from __future__ import division
from libtbx.test_utils import approx_equal
import iotbx.ncs
import iotbx.pdb

pdb_str_1 = """
CRYST1  399.000  399.000  399.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A  34     125.208 211.886 175.417  1.00  0.00           N
ATOM      2  CA  GLY A  34     125.035 211.123 174.168  1.00  0.00           C
ATOM      3  C   GLY A  34     126.386 210.806 173.507  1.00  0.00           C
ATOM      4  O   GLY A  34     127.304 211.628 173.503  1.00  0.00           O
TER
ATOM      5  N   GLY B  34     251.532 143.432 175.422  1.00  0.00           N
ATOM      6  CA  GLY B  34     252.120 143.948 174.173  1.00  0.00           C
ATOM      7  C   GLY B  34     251.212 144.998 173.512  1.00  0.00           C
ATOM      8  O   GLY B  34     249.986 144.872 173.510  1.00  0.00           O
TER
ATOM      9  N   GLY C  34     189.583 273.076 175.423  1.00  0.00           N
ATOM     10  CA  GLY C  34     188.804 273.006 174.173  1.00  0.00           C
ATOM     11  C   GLY C  34     188.920 271.622 173.510  1.00  0.00           C
ATOM     12  O   GLY C  34     189.986 271.004 173.508  1.00  0.00           O
TER
"""

pdb_str_2 = """
CRYST1  399.000  399.000  399.000  90.00  90.00  90.00 P 1
ATOM      1  O   HOH S   1     109.583 203.076 175.423  1.00  0.00           O
TER
ATOM      1  N   GLY A  34     125.208 211.886 175.417  1.00  0.00           N
ATOM      2  CA  GLY A  34     125.035 211.123 174.168  1.00  0.00           C
ATOM      3  C   GLY A  34     126.386 210.806 173.507  1.00  0.00           C
ATOM      4  O   GLY A  34     127.304 211.628 173.503  1.00  0.00           O
TER
ATOM      5  N   GLY B  34     251.532 143.432 175.422  1.00  0.00           N
ATOM      6  CA  GLY B  34     252.120 143.948 174.173  1.00  0.00           C
ATOM      7  C   GLY B  34     251.212 144.998 173.512  1.00  0.00           C
ATOM      8  O   GLY B  34     249.986 144.872 173.510  1.00  0.00           O
TER
ATOM      9  N   GLY C  34     189.583 273.076 175.423  1.00  0.00           N
ATOM     10  CA  GLY C  34     188.804 273.006 174.173  1.00  0.00           C
ATOM     11  C   GLY C  34     188.920 271.622 173.510  1.00  0.00           C
ATOM     12  O   GLY C  34     189.986 271.004 173.508  1.00  0.00           O
TER
ATOM      9  O   TYR D   4     189.583 273.076 175.423  1.00  0.00           O
TER
"""

pdb_str_3 = """
CRYST1  399.000  399.000  399.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A  34     125.208 211.886 175.417  1.00  0.00           N
ATOM      2  CA  GLY A  34     125.035 211.123 174.168  1.00  0.00           C
ATOM      3  C   GLY A  34     126.386 210.806 173.507  1.00  0.00           C
ATOM      4  O   GLY A  34     127.304 211.628 173.503  1.00  0.00           O
TER
ATOM      5  N   GLY B  34     251.532 143.432 175.422  1.00  0.00           N
ATOM      6  CA  GLY B  34     252.120 143.948 174.173  1.00  0.00           C
ATOM      7  C   GLY B  34     251.212 144.998 173.512  1.00  0.00           C
ATOM      8  O   GLY B  34     249.986 144.872 173.510  1.00  0.00           O
TER
ATOM      9  O   HOH C  34     189.583 273.076 175.423  1.00  0.00           O
ATOM     10  O   HOH C  34     188.804 273.006 174.173  1.00  0.00           O
ATOM     11  O   HOH C  34     188.920 271.622 173.510  1.00  0.00           O
ATOM     12  O   HOH C  34     189.986 271.004 173.508  1.00  0.00           O
TER
"""

pdb_str_4 = """
CRYST1  399.000  399.000  399.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A  34     125.208 211.886 175.417  1.00  0.00           N
ATOM      2  CA  GLY A  34     125.035 211.123 174.168  1.00  0.00           C
ATOM      3  C   GLY A  34     126.386 210.806 173.507  1.00  0.00           C
ATOM      4  O   GLY A  34     127.304 211.628 173.503  1.00  0.00           O
TER
ATOM      5  N   GLY B  34     251.532 143.432 175.422  1.00  0.00           N
ATOM      6  CA  GLY B  34     252.120 143.948 174.173  1.00  0.00           C
ATOM      7  C   GLY B  34     251.212 144.998 173.512  1.00  0.00           C
ATOM      8  O   GLY B  34     249.986 144.872 173.510  1.00  0.00           O
TER
ATOM      9  N   TYR C  34     189.583 273.076 175.423  1.00  0.00           N
ATOM     10  CA  TYR C  34     188.804 273.006 174.173  1.00  0.00           C
ATOM     11  C   TYR C  34     188.920 271.622 173.510  1.00  0.00           C
ATOM     12  O   TYR C  34     189.986 271.004 173.508  1.00  0.00           O
TER
"""

pdb_str_5 = """
ATOM      1  N   MET A   1     158.070 173.095 147.115  1.00 50.00           N
ATOM      2  CA  MET A   1     157.408 172.627 148.359  1.00 50.00           C
ATOM      3  CB  MET A   1     157.550 171.094 148.516  1.00 50.00           C
ATOM      4  CG  MET A   1     156.748 170.503 149.691  1.00 50.00           C
ATOM      5  SD  MET A   1     154.968 170.855 149.612  1.00 50.00           S
ATOM      6  CE  MET A   1     154.505 169.913 151.091  1.00 50.00           C
ATOM      7  C   MET A   1     157.958 173.331 149.563  1.00 50.00           C
ATOM      8  O   MET A   1     157.196 173.814 150.399  1.00 50.00           O
TER
ATOM      9  N   MET B   1     174.781 155.306 150.054  1.00 50.00           N
ATOM     10  CA  MET B   1     174.332 154.630 151.298  1.00 50.00           C
ATOM     11  CB  MET B   1     175.016 153.251 151.453  1.00 50.00           C
ATOM     12  CG  MET B   1     174.481 152.410 152.628  1.00 50.00           C
ATOM     13  SD  MET B   1     172.693 152.099 152.550  1.00 50.00           S
ATOM     14  CE  MET B   1     172.601 151.052 154.028  1.00 50.00           C
ATOM     15  C   MET B   1     174.594 155.484 152.502  1.00 50.00           C
ATOM     16  O   MET B   1     173.710 155.660 153.339  1.00 50.00           O
TER
ATOM     17  N   MET C   1     148.867 195.697 144.146  1.00 50.00           N
ATOM     18  CA  MET C   1     148.080 195.499 145.390  1.00 50.00           C
ATOM     19  CB  MET C   1     147.662 194.018 145.549  1.00 50.00           C
ATOM     20  CG  MET C   1     146.701 193.755 146.723  1.00 50.00           C
ATOM     21  SD  MET C   1     145.166 194.723 146.643  1.00 50.00           S
ATOM     22  CE  MET C   1     144.395 194.012 148.122  1.00 50.00           C
ATOM     23  C   MET C   1     148.846 195.960 146.594  1.00 50.00           C
ATOM     24  O   MET C   1     148.308 196.685 147.429  1.00 50.00           O
TER
ATOM    417  N   MET 1   1     274.499 237.478  69.907  1.00 50.00           N
ATOM    418  CA  MET 1   1     275.223 237.861  71.146  1.00 50.00           C
ATOM    419  CB  MET 1   1     275.281 239.400  71.298  1.00 50.00           C
ATOM    420  CG  MET 1   1     276.159 239.886  72.466  1.00 50.00           C
ATOM    421  SD  MET 1   1     277.878 239.307  72.379  1.00 50.00           S
ATOM    422  CE  MET 1   1     278.468 240.186  73.852  1.00 50.00           C
ATOM    423  C   MET 1   1     274.593 237.238  72.356  1.00 50.00           C
ATOM    424  O   MET 1   1     275.291 236.663  73.190  1.00 50.00           O
TER
END
"""

pdb_str_6 = """\
CRYST1  577.812  448.715  468.790  90.00  90.00  90.00 P 1
ATOM      1  CA  LYS A 151      10.766   9.333  12.905  1.00 44.22           C
ATOM      2  CA  LYS A 152      10.117   9.159  11.610  1.00 49.42           C
ATOM      3  CA  LYS A 153       9.099   8.000  11.562  1.00 46.15           C
ATOM      4  CA  LYS A 154       8.000   8.202  11.065  1.00 52.97           C
ATOM      5  CA  LYS A 155      11.146   9.065  10.474  1.00 41.68           C
ATOM      6  CA  LYS A 156      10.547   9.007   9.084  1.00 55.55           C
TER
ATOM      7  CA  LYS B 157      11.545   9.413   8.000  1.00 72.27           C
ATOM      8  CA  LYS B 158      12.277  10.718   8.343  1.00 75.78           C
ATOM      9  CA  LYS B 159      11.349  11.791   8.809  1.00 75.88           C
TER
ATOM      7  CA  LYS F 157       2.154   3.953  16.298  1.00 72.27           C
ATOM      8  CA  LYS F 158       2.014   3.732  17.811  1.00 75.78           C
ATOM      9  CA  LYS F 159       2.558   2.413  18.250  1.00 75.88           C
TER
ATOM      7  CA  LYS D 157       4.334  10.965  12.119  1.00 72.27           C
ATOM      8  CA  LYS D 158       4.057  11.980  13.238  1.00 75.78           C
ATOM      9  CA  LYS D 159       3.177  11.427  14.310  1.00 75.88           C
TER
ATOM      1  CA  LYS C 151       6.855   8.667  15.730  1.00 44.22           C
ATOM      2  CA  LYS C 152       5.891   8.459  14.655  1.00 49.42           C
ATOM      3  CA  LYS C 153       6.103   7.155  13.858  1.00 46.15           C
ATOM      4  CA  LYS C 154       5.138   6.438  13.633  1.00 52.97           C
ATOM      5  CA  LYS C 155       5.801   9.685  13.736  1.00 41.68           C
ATOM      6  CA  LYS C 156       4.731   9.594  12.667  1.00 55.55           C
TER
ATOM      1  CA  LYS E 151       6.987   4.106  17.432  1.00 44.22           C
ATOM      2  CA  LYS E 152       6.017   3.539  16.502  1.00 49.42           C
ATOM      3  CA  LYS E 153       6.497   3.492  15.036  1.00 46.15           C
ATOM      4  CA  LYS E 154       6.348   2.458  14.400  1.00 52.97           C
ATOM      5  CA  LYS E 155       4.647   4.221  16.634  1.00 41.68           C
ATOM      6  CA  LYS E 156       3.552   3.605  15.788  1.00 55.55           C
TER
"""

pdb_str_7 = """\
CRYST1  577.812  448.715  468.790  90.00  90.00  90.00 P 1
ATOM      1  CA  LYS A 151      10.766   9.333  12.905  1.00 44.22           C
ATOM      2  CA  LYS A 152      10.117   9.159  11.610  1.00 49.42           C
ATOM      3  CA  LYS A 153       9.099   8.000  11.562  1.00 46.15           C
ATOM      4  CA  LYS A 154       8.000   8.202  11.065  1.00 52.97           C
ATOM      5  CA  LYS A 155      11.146   9.065  10.474  1.00 41.68           C
TER
ATOM    222  CA  LEU X  40      94.618  -5.253  91.582  1.00 87.10           C
ATOM    223  CA  ARG X  41      62.395  51.344  80.786  1.00107.25           C
ATOM    224  CA  ARG X  42      62.395  41.344  80.786  1.00107.25           C
TER
ATOM      1  CA  THR D   1       8.111  11.080  10.645  1.00 20.00           C
ATOM      2  CA  THR D   2       8.000   9.722  10.125  1.00 20.00           C
ATOM      3  CA  THR D   3       8.075   8.694  11.249  1.00 20.00           C
ATOM      4  CA  THR D   4       8.890   8.818  12.163  1.00 20.00           C
TER
ATOM      1  CA  LYS B 151       6.855   8.667  15.730  1.00 44.22           C
ATOM      2  CA  LYS B 152       5.891   8.459  14.655  1.00 49.42           C
ATOM      3  CA  LYS B 153       6.103   7.155  13.858  1.00 46.15           C
ATOM      4  CA  LYS B 154       5.138   6.438  13.633  1.00 52.97           C
ATOM      5  CA  LYS B 155       5.801   9.685  13.736  1.00 41.68           C
TER
ATOM      1  CA  LYS C 151       6.987   4.106  17.432  1.00 44.22           C
ATOM      2  CA  LYS C 152       6.017   3.539  16.502  1.00 49.42           C
ATOM      3  CA  LYS C 153       6.497   3.492  15.036  1.00 46.15           C
ATOM      4  CA  LYS C 154       6.348   2.458  14.400  1.00 52.97           C
ATOM      5  CA  LYS C 155       4.647   4.221  16.634  1.00 41.68           C
TER
ATOM    222  CA  LEU Y  40     194.618   5.253  81.582  1.00 87.10           C
ATOM    223  CA  ARG Y  41     162.395  41.344  70.786  1.00107.25           C
ATOM    224  CA  ARG Y  42     162.395  31.344  70.786  1.00107.25           C
TER
ATOM      1  CA  THR E   1       8.111 -10.645  11.080  1.00 20.00           C
ATOM      2  CA  THR E   2       8.000 -10.125   9.722  1.00 20.00           C
ATOM      3  CA  THR E   3       8.075 -11.249   8.694  1.00 20.00           C
ATOM      4  CA  THR E   4       8.890 -12.163   8.818  1.00 20.00           C
TER
"""

pdb_str_8 = """\
ATOM      1  N   THR A   1      13.014  18.419   8.520  1.00 20.00           N
ATOM      2  CA  THR A   1      12.903  17.061   8.000  1.00 20.00           C
ATOM      3  C   THR A   1      12.978  16.033   9.124  1.00 20.00           C
ATOM      4  O   THR A   1      13.793  16.157  10.038  1.00 20.00           O
TER
ATOM      1  N   THR C   1      10.325   8.000  14.368  1.00 20.00           N
ATOM      2  CA  THR C   1      10.111   8.702  13.108  1.00 20.00           C
ATOM      3  C   THR C   1      11.313   9.570  12.750  1.00 20.00           C
ATOM      4  O   THR C   1      11.885  10.241  13.609  1.00 20.00           O
TER
ATOM      1  N   THR D   1      -0.430  15.458  18.495  1.00 20.00           N
ATOM      2  CA  THR D   1       0.086  15.020  17.204  1.00 20.00           C
ATOM      3  C   THR D   1       1.440  14.334  17.355  1.00 20.00           C
ATOM      4  O   THR D   1       2.297  14.791  18.111  1.00 20.00           O
TER
ATOM      1  N   THR E   1       1.895   8.366  19.752  1.00 20.00           N
ATOM      2  CA  THR E   1       1.682   9.068  18.491  1.00 20.00           C
ATOM      3  C   THR E   1       2.884   9.935  18.133  1.00 20.00           C
ATOM      4  O   THR E   1       3.455  10.606  18.993  1.00 20.00           O
TER
ATOM      1  N   THR F   1       8.346   7.308  15.936  1.00 20.00           N
ATOM      2  CA  THR F   1       7.054   7.796  15.467  1.00 20.00           C
ATOM      3  C   THR F   1       6.884   9.281  15.767  1.00 20.00           C
ATOM      4  O   THR F   1       7.237   9.751  16.849  1.00 20.00           O
TER
ATOM      1  N   THR G   1       0.609  -0.560  24.094  1.00 20.00           N
ATOM      2  CA  THR G   1       0.395   0.142  22.834  1.00 20.00           C
ATOM      3  C   THR G   1       1.598   1.009  22.476  1.00 20.00           C
ATOM      4  O   THR G   1       2.169   1.680  23.335  1.00 20.00           O
TER
ATOM      1  N   THR H   1       7.061  -1.617  20.279  1.00 20.00           N
ATOM      2  CA  THR H   1       5.768  -1.130  19.810  1.00 20.00           C
ATOM      3  C   THR H   1       5.599   0.356  20.109  1.00 20.00           C
ATOM      4  O   THR H   1       5.950   0.825  21.191  1.00 20.00           O
TER
ATOM      1  N   THR I   1       8.722   4.822  16.665  1.00 20.00           N
ATOM      2  CA  THR I   1       7.494   4.036  16.653  1.00 20.00           C
ATOM      3  C   THR I   1       6.628   4.350  17.868  1.00 20.00           C
ATOM      4  O   THR I   1       7.130   4.482  18.984  1.00 20.00           O
TER
ATOM      1  N   THR B   1       8.000  15.093  13.112  1.00 20.00           N
ATOM      2  CA  THR B   1       8.516  14.654  11.820  1.00 20.00           C
ATOM      3  C   THR B   1       9.870  13.968  11.972  1.00 20.00           C
ATOM      4  O   THR B   1      10.727  14.426  12.727  1.00 20.00           O
TER
"""

pdb_str_9 = """\
ATOM     45  N   PHEAa   6     221.693 146.930 114.416  1.00 50.00           N
ATOM     46  CA  PHEAa   6     220.871 148.020 114.886  1.00 50.00           C
ATOM     47  C   PHEAa   6     219.413 147.628 114.926  1.00 50.00           C
ATOM     48  O   PHEAa   6     218.730 147.905 115.908  1.00 50.00           O
ATOM     49  CB  PHEAa   6     221.058 149.265 113.976  1.00 50.00           C
ATOM     50  CG  PHEAa   6     220.338 150.481 114.498  1.00 50.00           C
ATOM     51  CD1 PHEAa   6     220.740 151.082 115.702  1.00 50.00           C
ATOM     52  CD2 PHEAa   6     219.240 151.016 113.801  1.00 50.00           C
ATOM     53  CE1 PHEAa   6     220.059 152.198 116.204  1.00 50.00           C
ATOM     54  CE2 PHEAa   6     218.557 152.133 114.299  1.00 50.00           C
ATOM     55  CZ  PHEAa   6     218.965 152.723 115.500  1.00 50.00           C
ATOM     56  N   ASNAa   7     218.926 146.940 113.868  1.00 50.00           N
ATOM     57  CA  ASNAa   7     217.565 146.462 113.741  1.00 50.00           C
ATOM     58  C   ASNAa   7     217.229 145.425 114.782  1.00 50.00           C
ATOM     59  O   ASNAa   7     216.113 145.412 115.283  1.00 50.00           O
ATOM     60  CB  ASNAa   7     217.235 145.875 112.347  1.00 50.00           C
ATOM     61  CG  ASNAa   7     217.277 146.970 111.271  1.00 50.00           C
ATOM     62  OD1 ASNAa   7     217.061 148.155 111.550  1.00 50.00           O
ATOM     63  ND2 ASNAa   7     217.551 146.541 110.002  1.00 50.00           N
ATOM     64  N   LEUAa   8     218.185 144.528 115.129  1.00 50.00           N
ATOM     65  CA  LEUAa   8     217.998 143.485 116.117  1.00 50.00           C
ATOM     66  C   LEUAa   8     217.846 144.051 117.505  1.00 50.00           C
ATOM     67  O   LEUAa   8     217.002 143.585 118.265  1.00 50.00           O
ATOM     68  CB  LEUAa   8     219.168 142.476 116.143  1.00 50.00           C
ATOM     69  CG  LEUAa   8     219.217 141.553 114.906  1.00 50.00           C
ATOM     70  CD1 LEUAa   8     220.600 140.889 114.770  1.00 50.00           C
ATOM     71  CD2 LEUAa   8     218.098 140.493 114.925  1.00 50.00           C
ATOM     72  N   LYSAa   9     218.665 145.069 117.866  1.00 50.00           N
ATOM     73  CA  LYSAa   9     218.573 145.769 119.133  1.00 50.00           C
ATOM     74  C   LYSAa   9     217.287 146.554 119.234  1.00 50.00           C
ATOM     75  O   LYSAa   9     216.682 146.606 120.301  1.00 50.00           O
ATOM     76  CB  LYSAa   9     219.751 146.744 119.364  1.00 50.00           C
ATOM     77  CG  LYSAa   9     221.113 146.057 119.558  1.00 50.00           C
ATOM     78  CD  LYSAa   9     221.257 145.305 120.891  1.00 50.00           C
ATOM     79  CE  LYSAa   9     222.646 144.673 121.065  1.00 50.00           C
ATOM     80  NZ  LYSAa   9     222.723 143.882 122.315  1.00 50.00           N
TER
ATOM   1244  N   PHEAb   1     305.367 162.705 105.239  1.00 50.00           N
ATOM   1245  CA  PHEAb   1     304.396 162.991 106.331  1.00 50.00           C
ATOM   1246  C   PHEAb   1     304.285 164.473 106.586  1.00 50.00           C
ATOM   1247  O   PHEAb   1     304.837 165.292 105.851  1.00 50.00           O
ATOM   1248  CB  PHEAb   1     304.743 162.176 107.627  1.00 50.00           C
ATOM   1249  CG  PHEAb   1     305.801 162.775 108.540  1.00 50.00           C
ATOM   1250  CD1 PHEAb   1     307.055 163.195 108.061  1.00 50.00           C
ATOM   1251  CD2 PHEAb   1     305.501 162.973 109.900  1.00 50.00           C
ATOM   1252  CE1 PHEAb   1     307.967 163.838 108.909  1.00 50.00           C
ATOM   1253  CE2 PHEAb   1     306.412 163.612 110.751  1.00 50.00           C
ATOM   1254  CZ  PHEAb   1     307.642 164.054 110.251  1.00 50.00           C
ATOM   1255  N   LYSAb   2     303.607 164.838 107.696  1.00 50.00           N
ATOM   1256  CA  LYSAb   2     303.542 166.190 108.180  1.00 50.00           C
ATOM   1257  C   LYSAb   2     303.642 166.054 109.667  1.00 50.00           C
ATOM   1258  O   LYSAb   2     302.934 165.255 110.277  1.00 50.00           O
ATOM   1259  CB  LYSAb   2     302.234 166.947 107.855  1.00 50.00           C
ATOM   1260  CG  LYSAb   2     302.035 167.201 106.352  1.00 50.00           C
ATOM   1261  CD  LYSAb   2     300.842 168.117 106.023  1.00 50.00           C
ATOM   1262  CE  LYSAb   2     301.067 169.582 106.424  1.00 50.00           C
ATOM   1263  NZ  LYSAb   2     299.910 170.421 106.031  1.00 50.00           N
ATOM   1264  N   ALAAb   3     304.533 166.867 110.282  1.00 50.00           N
ATOM   1265  CA  ALAAb   3     304.840 166.840 111.692  1.00 50.00           C
ATOM   1266  C   ALAAb   3     303.670 167.320 112.503  1.00 50.00           C
ATOM   1267  O   ALAAb   3     303.399 166.803 113.583  1.00 50.00           O
ATOM   1268  CB  ALAAb   3     306.055 167.711 112.038  1.00 50.00           C
ATOM   1269  N   GLUAb   4     302.961 168.344 111.975  1.00 50.00           N
ATOM   1270  CA  GLUAb   4     301.784 168.942 112.560  1.00 50.00           C
ATOM   1271  C   GLUAb   4     300.678 167.941 112.768  1.00 50.00           C
ATOM   1272  O   GLUAb   4     300.094 167.886 113.848  1.00 50.00           O
ATOM   1273  CB  GLUAb   4     301.233 170.108 111.709  1.00 50.00           C
ATOM   1274  CG  GLUAb   4     302.268 171.215 111.421  1.00 50.00           C
ATOM   1275  CD  GLUAb   4     302.808 171.795 112.727  1.00 50.00           C
ATOM   1276  OE1 GLUAb   4     304.040 171.677 112.965  1.00 50.00           O
ATOM   1277  OE2 GLUAb   4     301.994 172.363 113.504  1.00 50.00           O
TER
ATOM   2754  N   PHEAc   6     244.472 153.067 117.352  1.00 50.00           N
ATOM   2755  CA  PHEAc   6     243.314 153.789 117.823  1.00 50.00           C
ATOM   2756  C   PHEAc   6     242.094 152.900 117.864  1.00 50.00           C
ATOM   2757  O   PHEAc   6     241.358 152.912 118.847  1.00 50.00           O
ATOM   2758  CB  PHEAc   6     243.040 155.017 116.913  1.00 50.00           C
ATOM   2759  CG  PHEAc   6     241.932 155.894 117.436  1.00 50.00           C
ATOM   2760  CD1 PHEAc   6     242.091 156.600 118.639  1.00 50.00           C
ATOM   2761  CD2 PHEAc   6     240.715 155.999 116.739  1.00 50.00           C
ATOM   2762  CE1 PHEAc   6     241.055 157.396 119.141  1.00 50.00           C
ATOM   2763  CE2 PHEAc   6     239.676 156.795 117.237  1.00 50.00           C
ATOM   2764  CZ  PHEAc   6     239.846 157.493 118.439  1.00 50.00           C
ATOM   2765  N   ASNAc   7     241.886 152.082 116.806  1.00 50.00           N
ATOM   2766  CA  ASNAc   7     240.788 151.148 116.680  1.00 50.00           C
ATOM   2767  C   ASNAc   7     240.847 150.058 117.721  1.00 50.00           C
ATOM   2768  O   ASNAc   7     239.810 149.645 118.224  1.00 50.00           O
ATOM   2769  CB  ASNAc   7     240.690 150.481 115.286  1.00 50.00           C
ATOM   2770  CG  ASNAc   7     240.334 151.517 114.210  1.00 50.00           C
ATOM   2771  OD1 ASNAc   7     239.707 152.544 114.490  1.00 50.00           O
ATOM   2772  ND2 ASNAc   7     240.743 151.215 112.941  1.00 50.00           N
ATOM   2773  N   LEUAc   8     242.062 149.565 118.068  1.00 50.00           N
ATOM   2774  CA  LEUAc   8     242.264 148.525 119.056  1.00 50.00           C
ATOM   2775  C   LEUAc   8     241.919 148.999 120.443  1.00 50.00           C
ATOM   2776  O   LEUAc   8     241.299 148.260 121.204  1.00 50.00           O
ATOM   2777  CB  LEUAc   8     243.717 148.003 119.082  1.00 50.00           C
ATOM   2778  CG  LEUAc   8     244.094 147.160 117.844  1.00 50.00           C
ATOM   2779  CD1 LEUAc   8     245.623 147.037 117.707  1.00 50.00           C
ATOM   2780  CD2 LEUAc   8     243.430 145.768 117.864  1.00 50.00           C
ATOM   2781  N   LYSAc   9     242.317 150.242 120.804  1.00 50.00           N
ATOM   2782  CA  LYSAc   9     241.981 150.863 122.071  1.00 50.00           C
ATOM   2783  C   LYSAc   9     240.499 151.134 122.173  1.00 50.00           C
ATOM   2784  O   LYSAc   9     239.916 150.965 123.241  1.00 50.00           O
ATOM   2785  CB  LYSAc   9     242.729 152.197 122.302  1.00 50.00           C
ATOM   2786  CG  LYSAc   9     244.248 152.045 122.495  1.00 50.00           C
ATOM   2787  CD  LYSAc   9     244.653 151.394 123.828  1.00 50.00           C
ATOM   2788  CE  LYSAc   9     246.177 151.304 124.001  1.00 50.00           C
ATOM   2789  NZ  LYSAc   9     246.534 150.594 125.250  1.00 50.00           N
TER
ATOM   3953  N   PHEAd   1     316.882 197.854 108.123  1.00 50.00           N
ATOM   3954  CA  PHEAd   1     315.875 197.773 109.215  1.00 50.00           C
ATOM   3955  C   PHEAd   1     315.239 199.116 109.471  1.00 50.00           C
ATOM   3956  O   PHEAd   1     315.460 200.078 108.736  1.00 50.00           O
ATOM   3957  CB  PHEAd   1     316.493 197.136 110.512  1.00 50.00           C
ATOM   3958  CG  PHEAd   1     317.265 198.076 111.424  1.00 50.00           C
ATOM   3959  CD1 PHEAd   1     318.283 198.919 110.944  1.00 50.00           C
ATOM   3960  CD2 PHEAd   1     316.915 198.154 112.784  1.00 50.00           C
ATOM   3961  CE1 PHEAd   1     318.905 199.846 111.792  1.00 50.00           C
ATOM   3962  CE2 PHEAd   1     317.536 199.077 113.635  1.00 50.00           C
ATOM   3963  CZ  PHEAd   1     318.525 199.932 113.135  1.00 50.00           C
ATOM   3964  N   LYSAd   2     314.475 199.212 110.581  1.00 50.00           N
ATOM   3965  CA  LYSAd   2     313.929 200.450 111.066  1.00 50.00           C
ATOM   3966  C   LYSAd   2     314.073 200.360 112.553  1.00 50.00           C
ATOM   3967  O   LYSAd   2     313.699 199.360 113.163  1.00 50.00           O
ATOM   3968  CB  LYSAd   2     312.437 200.687 110.742  1.00 50.00           C
ATOM   3969  CG  LYSAd   2     312.159 200.852 109.239  1.00 50.00           C
ATOM   3970  CD  LYSAd   2     310.717 201.278 108.910  1.00 50.00           C
ATOM   3971  CE  LYSAd   2     310.400 202.726 109.311  1.00 50.00           C
ATOM   3972  NZ  LYSAd   2     309.019 203.094 108.919  1.00 50.00           N
ATOM   3973  N   ALAAd   3     314.612 201.440 113.166  1.00 50.00           N
ATOM   3974  CA  ALAAd   3     314.909 201.525 114.577  1.00 50.00           C
ATOM   3975  C   ALAAd   3     313.646 201.552 115.389  1.00 50.00           C
ATOM   3976  O   ALAAd   3     313.579 200.972 116.469  1.00 50.00           O
ATOM   3977  CB  ALAAd   3     315.731 202.774 114.922  1.00 50.00           C
ATOM   3978  N   GLUAd   4     312.615 202.253 114.862  1.00 50.00           N
ATOM   3979  CA  GLUAd   4     311.303 202.387 115.448  1.00 50.00           C
ATOM   3980  C   GLUAd   4     310.630 201.057 115.656  1.00 50.00           C
ATOM   3981  O   GLUAd   4     310.105 200.796 116.737  1.00 50.00           O
ATOM   3982  CB  GLUAd   4     310.369 203.279 114.597  1.00 50.00           C
ATOM   3983  CG  GLUAd   4     310.938 204.683 114.308  1.00 50.00           C
ATOM   3984  CD  GLUAd   4     311.233 205.418 115.614  1.00 50.00           C
ATOM   3985  OE1 GLUAd   4     312.425 205.752 115.850  1.00 50.00           O
ATOM   3986  OE2 GLUAd   4     310.271 205.657 116.392  1.00 50.00           O
TER
END
"""


def exercise_00(prefix="iotbx_ncs_exercise_00"):
  pdb_file_name = "%s.pdb"%prefix
  ncs_params_str = """
ncs_group {
  master_selection = chain A
  copy_selection = chain B
  copy_selection = chain C
}
  """
  def check_result(ncs_inp, test_i):
    if test_i == 0:
      l1, l2, l3 = [0,1,2,3], [4,5,6,7], [8,9,10,11]
    elif test_i == 1:
      l1, l2, l3 = [1,2,3,4], [5,6,7,8], [9,10,11,12]
    else: assert 0
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()
    assert len(ncs_groups) == 1
    ncs_group = ncs_groups[0]
    assert approx_equal(ncs_group.master_iselection, l1)
    assert len(ncs_group.copies) == 2
    assert approx_equal(ncs_group.copies[0].iselection, l2)
    assert approx_equal(ncs_group.copies[1].iselection, l3)
  for test_i, pdb_str in enumerate([pdb_str_1, pdb_str_2]):
    of = open(pdb_file_name, "w")
    print >> of, pdb_str
    of.close()
    pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
    if test_i == 0: # XXX Not implemented. Fix later.
      # using pdb_inp
      ncs_inp = iotbx.ncs.input(pdb_inp = pdb_inp)
      check_result(ncs_inp,test_i)
      # using file_name
      ncs_inp = iotbx.ncs.input(file_name = pdb_file_name)
      check_result(ncs_inp,test_i)
      # using pdb string
      ncs_inp = iotbx.ncs.input(pdb_string = pdb_str)
      check_result(ncs_inp,test_i)
    # using combination of pdb_inp and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_inp = pdb_inp,
      ncs_phil_string = ncs_params_str)
    check_result(ncs_inp,test_i)
    # using combination of pdb file name and Phil parameter string
    ncs_inp = iotbx.ncs.input(file_name = pdb_file_name,
      ncs_phil_string = ncs_params_str)
    check_result(ncs_inp,test_i)
    # using combination of pdb string and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_string = pdb_str,
      ncs_phil_string = ncs_params_str)
    check_result(ncs_inp,test_i)

def exercise_01(prefix="iotbx_ncs_exercise_01"):
  """
  Make sure provided selections take precedence and are correctly respected.
  """
  pdb_file_name = "%s.pdb"%prefix
  ncs_params_str = """
ncs_group {
  master_selection = chain C
  copy_selection = chain A
}
  """
  def check_result(ncs_inp, test_i):
    if test_i == 0:
      l1, l2 = [8,9,10,11], [0,1,2,3]
    elif test_i == 1:
      l1, l2 = [9,10,11,12], [1,2,3,4]
    else: assert 0
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()
    assert len(ncs_groups) == 1
    ncs_group = ncs_groups[0]
    assert approx_equal(ncs_group.master_iselection, l1)
    assert len(ncs_group.copies) == 1
    assert approx_equal(ncs_group.copies[0].iselection, l2)
  for test_i, pdb_str in enumerate([pdb_str_1, pdb_str_2]):
    of = open(pdb_file_name, "w")
    print >> of, pdb_str
    of.close()
    pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
    # using combination of pdb_inp and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_inp = pdb_inp,
      ncs_phil_string = ncs_params_str)
    check_result(ncs_inp,test_i)
    # using combination of pdb file name and Phil parameter string
    ncs_inp = iotbx.ncs.input(file_name = pdb_file_name,
      ncs_phil_string = ncs_params_str)
    check_result(ncs_inp,test_i)
    # using combination of pdb string and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_string = pdb_str,
      ncs_phil_string = ncs_params_str)
    check_result(ncs_inp,test_i)

def exercise_02(prefix="iotbx_ncs_exercise_02"):
  """
  This is expected to fail as requested chains cannot be matched.
  """
  pdb_file_name = "%s.pdb"%prefix
  ncs_params_str = """
ncs_group {
  master_selection = chain C
  copy_selection = chain A
}
  """
  for test_i, pdb_str in enumerate([pdb_str_3, pdb_str_4]):
    of = open(pdb_file_name, "w")
    print >> of, pdb_str
    of.close()
    pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
    # using combination of pdb_inp and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_inp = pdb_inp,
      ncs_phil_string = ncs_params_str)
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()
    # using combination of pdb file name and Phil parameter string
    ncs_inp = iotbx.ncs.input(file_name = pdb_file_name,
      ncs_phil_string = ncs_params_str)
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()
    # using combination of pdb string and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_string = pdb_str,
      ncs_phil_string = ncs_params_str)
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()

def exercise_03(prefix="iotbx_ncs_exercise_03"):
  """
  Expect one master and 3 copies.
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str_5)
  ncs_inp = iotbx.ncs.input(pdb_inp = pdb_inp)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  asc = ncs_inp.hierarchy.atom_selection_cache()
  sel_master = asc.selection(string = "chain 1")
  sel_copy_1 = asc.selection(string = "chain A")
  sel_copy_2 = asc.selection(string = "chain B")
  sel_copy_3 = asc.selection(string = "chain C")
  assert len(ncs_groups)==1
  ng = ncs_groups[0]
  # chains are sorted by name (numbers first)
  assert approx_equal(sel_master.iselection(), ng.master_iselection)
  assert approx_equal(sel_copy_1.iselection(), ng.copies[0].iselection)
  assert approx_equal(sel_copy_2.iselection(), ng.copies[1].iselection)
  assert approx_equal(sel_copy_3.iselection(), ng.copies[2].iselection)

def exercise_04(prefix="iotbx_ncs_exercise_04"):
  """
  Testing one ncs group and master ncs with two chains of different length
  and chains are not in order
  """
  ncs_inp = iotbx.ncs.input(pdb_string=pdb_str_6)
  t = ncs_inp.ncs_to_asu_selection
  assert t.keys() == ['chain A or chain B']
  assert t.values() == [['chain C or chain D', 'chain E or chain F']]

def exercise_05(prefix="iotbx_ncs_exercise_05"):
  """
  Make sure that phil selection overrides ncs grouping
  """
  phil_str = """\
ncs_group {
  master_selection = chain A
  copy_selection = chain C
}
ncs_group {
  master_selection = chain B
  copy_selection = chain D
}
"""
  ncs_inp = iotbx.ncs.input(
    pdb_string=pdb_str_6,
    ncs_phil_string=phil_str)
  expected = {'chain A': ['chain C'], 'chain B': ['chain D']}
  assert ncs_inp.ncs_to_asu_selection.keys(), expected.keys()
  assert ncs_inp.ncs_to_asu_selection.values(), expected.values()

def exercise_06(prefix="iotbx_ncs_exercise_06"):
  """
  Two groups, different number of chain in each group.
  Two chains that are NOT ncs related
  """
  ncs_inp = iotbx.ncs.input(pdb_string=pdb_str_7)
  t = ncs_inp.ncs_to_asu_selection
  assert t.keys() == ['chain D', 'chain A']
  assert t.values() == [['chain E'], ['chain B', 'chain C']]
  assert ncs_inp.ncs_group_map[1][0] == {'chain A'}
  assert ncs_inp.ncs_group_map[2][0] == {'chain D'}

def exercise_07(prefix="iotbx_ncs_exercise_07"):
  """
  Test that minimal number of chains in master ncs are selected (not the
  minimal number of transformations)
  """
  ncs_inp = iotbx.ncs.input(pdb_string=pdb_str_8)
  t = ncs_inp.ncs_to_asu_selection
  assert t.keys() == ['chain A']
  assert t.values() == \
    [['chain B', 'chain C', 'chain D', 'chain E',
      'chain F', 'chain G', 'chain H', 'chain I']]

def exercise_08(prefix="iotbx_ncs_exercise_08"):
  """
  Test that minimal number of transformations in master ncs are selected
  (not the minimal number of chains)
  """
  ncs_inp = iotbx.ncs.input(
    pdb_string=pdb_str_8,
    use_minimal_master_ncs=False)
  t = ncs_inp.ncs_to_asu_selection
  assert t.keys()==['chain A or chain B or chain C']
  assert t.values()==\
    [['chain D or chain E or chain F',
      'chain G or chain H or chain I']]

def exercise_09(prefix="iotbx_ncs_exercise_09"):
  """ ??? """
  from mmtbx.ncs import ncs_search
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str_8)
  ph = pdb_inp.construct_hierarchy()
  #
  chain_match_list = ncs_search.search_ncs_relations(
    ph=ph,min_contig_length=0,min_percent=0)
  # make sure that all possible chains compared
  model  = ph.models()[0]
  chain_ids = list({x.id for x in model.chains()})
  chain_ids = sorted(chain_ids)
  n_chains = len(chain_ids)
  assert n_chains==9
  assert len(chain_match_list)==8
  #
  match_dict = ncs_search.clean_chain_matching(chain_match_list,ph)
  chains_info = ncs_search.get_chains_info(ph)
  # Test minimal master NCS
  transform_to_group,match_dict = ncs_search.minimal_master_ncs_grouping(
  match_dict=match_dict)
  group_dict = ncs_search.build_group_dict(
    transform_to_group,match_dict,chains_info)
  assert len(group_dict)==1
  gr_obj = group_dict[('A',)]
  assert len(gr_obj.transforms)==len(gr_obj.copies)
  assert len(gr_obj.iselections)==len(gr_obj.copies)
  expected = [['A'], ['B'], ['C'], ['D'], ['E'], ['F'], ['G'], ['H'], ['I']]
  assert gr_obj.copies == expected
  tr = gr_obj.transforms[0]
  assert tr.r.is_r3_identity_matrix()
  assert tr.t.is_col_zero()
  tr = gr_obj.transforms[1]
  assert not tr.r.is_r3_identity_matrix()
  assert not tr.t.is_col_zero()

def exercise_10(prefix="iotbx_ncs_exercise_10"):
  """ Test minimal NCS operators """
  from mmtbx.ncs import ncs_search
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str_8)
  ph = pdb_inp.construct_hierarchy()
  #
  chain_match_list = ncs_search.search_ncs_relations(
    ph=ph,min_contig_length=0,
    min_percent=0,
    use_minimal_master_ncs=False)
  #
  match_dict = ncs_search.clean_chain_matching(chain_match_list,ph)
  chains_info = ncs_search.get_chains_info(ph)
  #
  transform_to_group,match_dict = ncs_search.minimal_ncs_operators_grouping(
  match_dict=match_dict)
  group_dict = ncs_search.build_group_dict(
    transform_to_group,match_dict,chains_info)
  assert len(group_dict)==1
  gr_obj = group_dict[('A', 'B', 'C')]
  assert len(gr_obj.transforms)==len(gr_obj.copies)
  assert len(gr_obj.iselections)==len(gr_obj.copies)
  expected = [['A', 'B', 'C'], ['D', 'E', 'F'], ['G', 'H', 'I']]
  assert gr_obj.copies==expected
  tr = gr_obj.transforms[0]
  assert tr.r.is_r3_identity_matrix()
  assert tr.t.is_col_zero()
  tr = gr_obj.transforms[1]
  assert not tr.r.is_r3_identity_matrix()
  assert not tr.t.is_col_zero()

def exercise_11(prefix="iotbx_ncs_exercise_11"):
  """
  Make sure user-provided NCS groups are preserved. Also make sure two-letter
  chain ID are handled correctly
  """
  phil_str="""
ncs_group {
  master_selection = chain 'Aa' and (resseq 6:9 )
  copy_selection = chain 'Ac' and (resseq 6:9 )
}
ncs_group {
  master_selection = chain 'Ab' and (resseq 1:4 )
  copy_selection = chain 'Ad' and (resseq 1:4 )
}
"""
  ncs_inp = iotbx.ncs.input(pdb_string = pdb_str_9)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  assert len(ncs_groups)==1
  #
  ncs_inp = iotbx.ncs.input(pdb_string = pdb_str_9, ncs_phil_string = phil_str)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  assert len(ncs_groups)==2

if (__name__ == "__main__"):
  exercise_00()
  exercise_01()
  exercise_02()
  exercise_03()
  exercise_04()
  exercise_05()
  exercise_06()
  exercise_07()
  exercise_08()
  exercise_09()
  exercise_10()
  exercise_11()
